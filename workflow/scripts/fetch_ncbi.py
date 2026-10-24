import os
import io
import re
import sys
import yaml
import random
import pandas as pd
from Bio import Entrez, SeqIO

def clean_accession(raw_id):
    """Extracts clean base accession from NCBI or ENA."""
    cleaned = str(raw_id).strip().strip("|")
    if "|" in cleaned:
        cleaned = cleaned.split("|")[-1]
    match = re.search(r"([A-Z]{1,6}_?[0-9]{4,})", cleaned)
    return match.group(1) if match else cleaned.split(".")[0]

def get_nucleotide_ids_for_assembly(assembly_acc):
    """Finds linked nucleotide UIDs for a GCA/GCF assembly without duplicating RefSeq+GenBank or grabbing WGS masters."""
    search_handle = Entrez.esearch(db="assembly", term=f"{assembly_acc}[Assembly Accession]")
    search_record = Entrez.read(search_handle)
    if not search_record["IdList"]:
        # Fallback to unscoped term if needed
        search_handle = Entrez.esearch(db="assembly", term=assembly_acc)
        search_record = Entrez.read(search_handle)
        if not search_record["IdList"]:
            return None

    asm_id = search_record["IdList"][0]
    link_handle = Entrez.elink(dbfrom="assembly", db="nucleotide", id=asm_id)
    link_record = Entrez.read(link_handle)

    try:
        link_sets = {
            ls["LinkName"]: [link["Id"] for link in ls["Link"]]
            for ls in link_record[0].get("LinkSetDb", [])
        }
        # Prefer RefSeq for GCF_ and INSDC (GenBank) for GCA_; never use wgsmaster
        if assembly_acc.upper().startswith("GCF_"):
             preferred_order = ["assembly_nuccore_refseq", "assembly_nuccore_insdc"]
        else:
             preferred_order = ["assembly_nuccore_insdc", "assembly_nuccore_refseq"]

        for name in preferred_order:
            if name in link_sets and link_sets[name]:
                return link_sets[name]

        # Fallback to any non-wgsmaster linkset
        for name, ids in link_sets.items():
            if "wgsmaster" not in name and ids:
                return ids
    except (IndexError, KeyError):
        pass
    return None

def main():
    with open("config/config.yaml", "r") as f:
        config = yaml.safe_load(f)

    os.makedirs("results/raw", exist_ok=True)
    Entrez.email = config.get("ncbi_email") or "test@example.com"

    ref_acc = str(config["reference_accession"]).strip()
    base_ref = clean_accession(ref_acc)
    print(f"🧬 Fetching Reference Genome & GenBank Annotations: {ref_acc}")

    fetch_id = ref_acc
    if ref_acc.upper().startswith(("GCF_", "GCA_")):
        print(f"🔗 Assembly detected. Resolving {ref_acc} to nucleotide sequences...")
        nuc_uids = get_nucleotide_ids_for_assembly(ref_acc)
        if not nuc_uids:
            print(f"❌ Failed to find nucleotide sequences for assembly {ref_acc}", file=sys.stderr)
            sys.exit(1)
        fetch_id = ",".join(nuc_uids)

    try:
        with Entrez.efetch(db="nucleotide", id=fetch_id, rettype="gbwithparts", retmode="text") as handle:
            gb_text = handle.read()

        with open("results/raw/reference.gb", "w") as gb_out:
            gb_out.write(gb_text)

        ref_records = list(SeqIO.parse(io.StringIO(gb_text), "genbank"))
        if not ref_records:
            raise ValueError("No GenBank records found in reference handle.")

        original_ref_accs = []
        for rec in ref_records:
            cid = clean_accession(rec.id)
            original_ref_accs.append(cid)
            rec.id = cid
            rec.name = cid
            rec.description = ""

        SeqIO.write(ref_records, "results/raw/reference.fasta", "fasta")
        max_ref_len = max(len(r.seq) for r in ref_records)
        print(f"✅ Loaded {len(ref_records)} reference segment(s)/contig(s) (largest: {max_ref_len:,} bp).")

    except Exception as e:
        print(f"❌ Failed to fetch NCBI Reference: {e}", file=sys.stderr)
        sys.exit(1)

    # 2. Get a pool of wild isolates
    search_mode = config.get("search_mode", "By Organism & Keyword")
    ncbi_query = config.get("search_term", '"Respiratory syncytial virus"[Organism]')
    min_len = int(config.get("min_length", 10000))
    max_len = max(50000, int(max_ref_len * 1.5))
    max_n_percent = float(config.get("max_ambiguous_n_percent", 1.0))
    target_quota = int(config.get("max_sequences", 20))

    if search_mode == "By Accession List":
        raw_accs = [x.strip() for x in re.split(r"[,\s]+", str(config.get("accession_list", ""))) if x.strip()]
        id_list = []
        for acc in raw_accs:
            if acc.upper().startswith(("GCF_", "GCA_")):
                asm_uids = get_nucleotide_ids_for_assembly(acc)
                if asm_uids:
                    id_list.extend(asm_uids)
            else:
                with Entrez.esearch(db="nucleotide", term=f"{acc}[Accession]", retmax=5) as h:
                    res = Entrez.read(h)
                    id_list.extend(res.get("IdList", []))
        target_quota = max(target_quota, len(id_list))
        print(f"🔍 Resolved {len(id_list)} nucleotide record(s) from Accession List.")
    else:
        exclude_accs = original_ref_accs[:10]
        exclude_str = " ".join([f"NOT {acc}[Accession]" for acc in exclude_accs])
        full_query = f"({ncbi_query}) AND {min_len}:{max_len}[SLEN] {exclude_str}".strip()
        print(f"🔍 Executing NCBI Query: {full_query}")

        with Entrez.esearch(db="nucleotide", term=full_query, retmax=5000) as handle:
            record = Entrez.read(handle)
        id_list = record.get("IdList", [])
        random.shuffle(id_list)

    print(f"📦 Found {len(id_list)} candidate sequences in NCBI.")
    if not id_list:
        print("❌ NCBI returned no sequences for query.", file=sys.stderr)
        sys.exit(1)

    valid_records = []
    seen_accessions = set(original_ref_accs + [base_ref])
    chunk_size = 30

    print(f"🧬 Actively hunting for {target_quota} pristine sequences...")

    for i in range(0, len(id_list), chunk_size):
        chunk_uids = id_list[i : i + chunk_size]
        try:
            with Entrez.efetch(db="nucleotide", id=",".join(chunk_uids), rettype="fasta", retmode="text") as fetch_handle:
                records = list(SeqIO.parse(fetch_handle, "fasta"))
        except Exception:
            continue

        for rec in records:
            clean_id = clean_accession(rec.id)
            seq_str = str(rec.seq).upper()

            if len(seq_str) == 0:
                continue

            n_count = seq_str.count("N")
            n_percent = (n_count / len(seq_str)) * 100

            if clean_id not in seen_accessions and len(seq_str) >= min_len and n_percent <= max_n_percent:
                seen_accessions.add(clean_id)
                rec.id = clean_id
                rec.name = clean_id
                # Preserve original title in a custom attribute before clearing description for FASTA
                rec.annotations["orig_title"] = rec.description
                rec.description = ""
                valid_records.append(rec)

            if len(valid_records) >= target_quota:
                break

        if len(valid_records) >= target_quota:
            break

    if not valid_records:
        print("❌ No sequences passed the QC filters during fetching.", file=sys.stderr)
        sys.exit(1)

    SeqIO.write(valid_records, "results/raw/isolates.fasta", "fasta")

    # Fetch Metadata by exact accession so skipped records never shift UIDs
    valid_accs = [r.id for r in valid_records]
    print(f"📊 Fetching metadata for the {len(valid_accs)} successful sequences...")
    meta_by_acc = {}
    try:
        with Entrez.esummary(db="nucleotide", id=",".join(valid_accs)) as summary_handle:
            summaries = Entrez.read(summary_handle)
        for s in summaries:
            raw_acc = s.get("Caption") or s.get("AccessionVersion") or ""
            clean_acc = clean_accession(raw_acc)
            meta_by_acc[clean_acc] = {
                "Accession": clean_acc,
                "Title": s.get("Title", ""),
                "Length_bp": int(s.get("Length", 0)),
                "UpdateDate": s.get("UpdateDate", ""),
            }
    except Exception as e:
        print(f"⚠️ Warning: esummary fallback triggered ({e})")

    meta_data = []
    for r in valid_records:
        if r.id in meta_by_acc:
            meta_data.append(meta_by_acc[r.id])
        else:
            meta_data.append({
                "Accession": r.id,
                "Title": r.annotations.get("orig_title", r.id),
                "Length_bp": len(r.seq),
                "UpdateDate": "",
            })

    pd.DataFrame(meta_data).to_csv("results/raw/metadata.csv", index=False)
    print(f"✅ Success! Saved {len(valid_records)} pristine NCBI sequences and metadata.")

if __name__ == "__main__":
    main()