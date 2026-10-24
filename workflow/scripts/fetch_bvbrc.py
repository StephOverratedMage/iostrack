import os
import io
import re
import sys
import yaml
import random
import requests
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def clean_accession(raw_id):
    cleaned = str(raw_id).strip().strip("|")
    if "|" in cleaned:
        cleaned = cleaned.split("|")[-1]
    match = re.search(r"([A-Z]{1,6}_?[0-9]{4,})", cleaned)
    return match.group(1) if match else cleaned.split(".")[0]

def get_nucleotide_ids_for_assembly(assembly_acc):
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

    # 1. Fetch Reference Genome from NCBI (required for rich GenBank CDS features)
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

    # 2. Build BV-BRC Query
    bvbrc_query = config.get("bvbrc_query", "eq(taxon_id,11250)")
    target_quota = int(config.get("max_sequences", 20))
    min_len = int(config.get("min_length", 10000))
    max_len = max(50000, int(max_ref_len * 1.5))
    max_n_percent = float(config.get("max_ambiguous_n_percent", 1.0))
    search_mode = config.get("search_mode", "By Organism & Keyword")

    if search_mode == "By Accession List":
        full_query = f"{bvbrc_query}&select(accession,description,sequence,length)"
    else:
        full_query = f"{bvbrc_query}&ge(length,{min_len})&le(length,{max_len})&select(accession,description,sequence,length)&limit({max(target_quota * 5, 200)})"

    print(f"🌍 Querying BV-BRC API for candidate pool: {full_query}")

    url = "https://www.bv-brc.org/api/genome_sequence/"
    headers = {
        "Accept": "application/json",
        "Content-Type": "application/x-www-form-urlencoded"
    }

    response = requests.post(url, data=full_query, headers=headers, timeout=60)
    if response.status_code != 200:
        print(f"❌ BV-BRC API Error {response.status_code}: {response.text}", file=sys.stderr)
        sys.exit(1)

    records_data = response.json()
    if not records_data:
        print("❌ BV-BRC returned HTTP 200 but NO sequences were found.", file=sys.stderr)
        sys.exit(1)
        
    random.shuffle(records_data)
    print(f"🧬 Actively hunting for {target_quota} pristine sequences from {len(records_data)} BV-BRC candidates...")

    valid_records = []
    meta_data = []
    seen_accessions = set(original_ref_accs + [base_ref])

    for item in records_data:
        raw_acc = item.get("accession", "")
        clean_id = clean_accession(raw_acc)
        seq_str = str(item.get("sequence", "")).upper()
        
        if not seq_str:
            continue
            
        n_count = seq_str.count("N")
        n_percent = (n_count / len(seq_str)) * 100

        if clean_id not in seen_accessions and len(seq_str) >= min_len and n_percent <= max_n_percent:
            seen_accessions.add(clean_id)
            
            # Create SeqRecord
            rec = SeqRecord(
                Seq(seq_str),
                id=clean_id,
                name=clean_id,
                description=""
            )
            valid_records.append(rec)
            
            # Store metadata
            meta_data.append({
                "Accession": clean_id,
                "Title": item.get("description", raw_acc),
                "Length_bp": int(item.get("length", len(seq_str))),
                "UpdateDate": ""  # BV-BRC sequence API omits upload dates
            })
            
        if len(valid_records) >= target_quota:
            break

    if not valid_records:
        print("❌ No BV-BRC sequences passed the QC filters.", file=sys.stderr)
        sys.exit(1)

    SeqIO.write(valid_records, "results/raw/isolates.fasta", "fasta")
    pd.DataFrame(meta_data).to_csv("results/raw/metadata.csv", index=False)
    print(f"✅ Success! Saved {len(valid_records)} pristine BV-BRC sequences and metadata.")

if __name__ == "__main__":
    main()