import os
import io
import re
import sys
import yaml
import requests
import pandas as pd
from Bio import Entrez, SeqIO


def clean_accession(raw_id):
    """Extracts clean base accession from NCBI or ENA."""
    cleaned = str(raw_id).strip().strip("|")
    if "|" in cleaned:
        cleaned = cleaned.split("|")[-1]
    match = re.search(r"([A-Z]{1,6}_?[0-9]{4,})", cleaned)
    return match.group(1) if match else cleaned.split(".")[0]


def load_local_fasta_isolates(config, seen_accessions, min_len, max_n_percent):
    """Loads user-uploaded local FASTA sequences and applies QC filters."""
    local_path = config.get("local_fasta_path", "")
    local_records = []
    local_meta = []

    if not local_path or not os.path.exists(local_path):
        return local_records, local_meta

    print(f"📂 Reading uploaded local FASTA file: {local_path}")
    for rec in SeqIO.parse(local_path, "fasta"):
        raw_header = rec.description or rec.id
        sanitized = re.sub(r"[^A-Za-z0-9_.-]", "_", str(rec.id).strip())
        clean_id = clean_accession(sanitized)
        seq_str = str(rec.seq).upper()
        if not seq_str:
            continue

        n_percent = (seq_str.count("N") / len(seq_str)) * 100.0
        if clean_id not in seen_accessions and len(seq_str) >= min_len and n_percent <= max_n_percent:
            seen_accessions.add(clean_id)
            rec.id = clean_id
            rec.name = clean_id
            rec.description = ""
            local_records.append(rec)
            local_meta.append({
                "Accession": clean_id,
                "Title": f"[Local FASTA] {raw_header}",
                "Length_bp": len(seq_str),
                "UpdateDate": "Local",
            })

    print(f"✅ Accepted {len(local_records)} local FASTA sequence(s) after QC.")
    return local_records, local_meta


def get_nucleotide_ids_for_assembly(assembly_acc):
    """Finds linked nucleotide UIDs for a GCA/GCF assembly without duplicating RefSeq+GenBank or grabbing WGS masters."""
    search_handle = Entrez.esearch(db="assembly", term=f"{assembly_acc}[Assembly Accession]")
    search_record = Entrez.read(search_handle)
    if not search_record["IdList"]:
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
        if assembly_acc.upper().startswith("GCF_"):
            preferred_order = ["assembly_nuccore_refseq", "assembly_nuccore_insdc"]
        else:
            preferred_order = ["assembly_nuccore_insdc", "assembly_nuccore_refseq"]

        for name in preferred_order:
            if name in link_sets and link_sets[name]:
                return link_sets[name]

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

    # 2. Build ENA Query & Fetch Candidate Metadata Pool
    search_mode = config.get("search_mode", "By Organism & Keyword")
    ena_query = config.get("ena_query", "tax_tree(11250)")
    target_quota = int(config.get("max_sequences", 20))
    min_len = int(config.get("min_length", 10000))
    max_len = max(50000, int(max_ref_len * 1.5))
    max_n_percent = float(config.get("max_ambiguous_n_percent", 1.0))

    seen_accessions = set(original_ref_accs + [base_ref])
    local_records, local_meta = load_local_fasta_isolates(
        config, seen_accessions, min_len, max_n_percent
    )

    if search_mode == "By Accession List":
        raw_accs = [x.strip() for x in re.split(r"[,\s]+", str(config.get("accession_list", ""))) if x.strip()]
        target_quota = max(target_quota, len(raw_accs))
        full_query = ena_query
    else:
        full_query = f"{ena_query} AND base_count>={min_len} AND base_count<={max_len}"
    print(f"🌍 Querying ENA Portal API for candidate pool: {full_query}")

    portal_url = "https://www.ebi.ac.uk/ena/portal/api/search"
    meta_params = {
        "result": "sequence",
        "query": full_query,
        "limit": max(target_quota * 5, 200),
        "format": "tsv",
        "fields": "accession,description,base_count,last_updated",
    }

    meta_response = requests.get(portal_url, params=meta_params, timeout=30)
    if meta_response.status_code != 200:
        print(f"❌ ENA Metadata Error {meta_response.status_code}: {meta_response.text}", file=sys.stderr)
        sys.exit(1)

    meta_text = meta_response.text.strip()
    if not meta_text and not local_records:
        print(f"❌ ENA returned HTTP 200 but NO sequences were found for query: {full_query}", file=sys.stderr)
        sys.exit(1)

    valid_records = list(local_records)
    valid_ids = set()
    total_goal = target_quota + len(local_records)

    if meta_text:
        meta_df = pd.read_csv(io.StringIO(meta_text), sep="\t")
        if "accession" in meta_df.columns and not meta_df.empty:
            meta_df["accession"] = meta_df["accession"].astype(str).apply(clean_accession)
            meta_df = meta_df[~meta_df["accession"].isin(seen_accessions)].drop_duplicates(subset=["accession"])
            meta_df = meta_df.sample(frac=1.0).reset_index(drop=True)
            candidate_ids = meta_df["accession"].tolist()

            print(f"🧬 Actively hunting for {target_quota} pristine ENA sequences...")
            chunk_size = 30

            for i in range(0, len(candidate_ids), chunk_size):
                chunk = candidate_ids[i : i + chunk_size]
                browser_url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{','.join(chunk)}"
                fasta_response = requests.get(browser_url, timeout=30)

                if fasta_response.status_code != 200:
                    continue

                records = list(SeqIO.parse(io.StringIO(fasta_response.text), "fasta"))
                for record in records:
                    clean_id = clean_accession(record.id)
                    seq_str = str(record.seq).upper()

                    if len(seq_str) == 0:
                        continue

                    n_count = seq_str.count("N")
                    n_percent = (n_count / len(seq_str)) * 100

                    if clean_id not in seen_accessions and len(seq_str) >= min_len and n_percent <= max_n_percent:
                        seen_accessions.add(clean_id)
                        valid_ids.add(clean_id)
                        record.id = clean_id
                        record.name = clean_id
                        record.description = ""
                        valid_records.append(record)

                    if len(valid_records) >= total_goal:
                        break
                if len(valid_records) >= total_goal:
                    break

            final_meta_df = meta_df[meta_df["accession"].isin(valid_ids)].copy()
            final_meta_df = final_meta_df.rename(
                columns={"accession": "Accession", "description": "Title", "base_count": "Length_bp", "last_updated": "UpdateDate"}
            )
        else:
            final_meta_df = pd.DataFrame(columns=["Accession", "Title", "Length_bp", "UpdateDate"])
    else:
        final_meta_df = pd.DataFrame(columns=["Accession", "Title", "Length_bp", "UpdateDate"])

    if not valid_records:
        print("❌ No ENA or local sequences passed the QC filters.", file=sys.stderr)
        sys.exit(1)

    SeqIO.write(valid_records, "results/raw/isolates.fasta", "fasta")

    if local_meta:
        final_meta_df = pd.concat([pd.DataFrame(local_meta), final_meta_df], ignore_index=True)

    final_meta_df.to_csv("results/raw/metadata.csv", index=False)
    print(f"✅ Success! Saved {len(valid_records)} pristine sequence(s) and metadata.")


if __name__ == "__main__":
    main()