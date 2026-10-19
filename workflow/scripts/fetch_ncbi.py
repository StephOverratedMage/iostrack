import os
import io
import re
import sys
import yaml
import random
import pandas as pd
from Bio import Entrez, SeqIO

def clean_accession(raw_id):
    """Extracts clean base accession from NCBI or ENA (e.g. ENA|OR287049|OR287049.1 -> OR287049)."""
    cleaned = str(raw_id).strip().strip("|")
    if "|" in cleaned:
        cleaned = cleaned.split("|")[-1]
    match = re.search(r'([A-Z]{1,6}_?[0-9]{4,})', cleaned)
    return match.group(1) if match else cleaned.split(".")[0]

def main():
    with open("config/config.yaml", "r") as f:
        config = yaml.safe_load(f)
        
    os.makedirs("results/raw", exist_ok=True)
    Entrez.email = config.get("ncbi_email", "test@example.com")
    
    # 1. Fetch Reference Genome (GenBank + FASTA in one call)
    ref_acc = config["reference_accession"]
    base_ref = clean_accession(ref_acc)
    Entrez.email = config.get("ncbi_email", "test@example.com")
    
    print(f"🧬 Fetching Reference Genome & GenBank Annotations: {ref_acc}")
    try:
        with Entrez.efetch(db="nucleotide", id=ref_acc, rettype="gbwithparts", retmode="text") as handle:
            gb_text = handle.read()
            
        with open("results/raw/reference.gb", "w") as gb_out:
            gb_out.write(gb_text)
            
        ref_record = SeqIO.read(io.StringIO(gb_text), "genbank")
        ref_record.id = base_ref
        ref_record.name = base_ref
        ref_record.description = ""
        SeqIO.write(ref_record, "results/raw/reference.fasta", "fasta")
    except Exception as e:
        print(f"❌ Failed to fetch NCBI Reference: {e}", file=sys.stderr)
        sys.exit(1)
        
    # 2. Get a pool of wild isolates
    ncbi_query = config.get("search_term", '"Respiratory Syncytial Virus"[Organism]')
    min_len = int(config.get("min_length", 10000))
    max_n_percent = float(config.get("max_ambiguous_n_percent", 1.0))
    target_quota = int(config.get("max_sequences", 20))
    
    # Filter by length in Entrez and exclude the reference accession
    full_query = f"({ncbi_query}) AND {min_len}:50000[SLEN] NOT {base_ref}[Accession]"
    print(f"🔍 Executing NCBI Query: {full_query}")
    
    with Entrez.esearch(db="nucleotide", term=full_query, retmax=5000) as handle:
        record = Entrez.read(handle)
        
    id_list = record.get("IdList", [])
    print(f"📦 Found {len(id_list)} candidate sequences in NCBI.")
    if not id_list:
        print("❌ NCBI returned no sequences for query.", file=sys.stderr)
        sys.exit(1)
        
    random.shuffle(id_list)
    
    # 3. ACTIVE VALIDATION LOOP: Fetch in chunks and retain numeric UIDs
    valid_records = []
    valid_uids = []
    seen_accessions = {base_ref}
    chunk_size = 30
    
    print(f"🧬 Actively hunting for {target_quota} pristine sequences...")
    
    for i in range(0, len(id_list), chunk_size):
        chunk_uids = id_list[i:i+chunk_size]
        try:
            with Entrez.efetch(db="nucleotide", id=",".join(chunk_uids), rettype="fasta", retmode="text") as fetch_handle:
                records = list(SeqIO.parse(fetch_handle, "fasta"))
        except Exception as e:
            print(f"⚠️ Warning: Chunk fetch failed ({e}), skipping chunk...")
            continue
            
        # Pair each numeric UID with its downloaded FASTA record
        for uid, rec in zip(chunk_uids, records):
            clean_id = clean_accession(rec.id)
            seq_str = str(rec.seq).upper()
            
            if len(seq_str) == 0:
                continue
                
            n_count = seq_str.count('N')
            n_percent = (n_count / len(seq_str)) * 100
            
            if clean_id not in seen_accessions and len(seq_str) >= min_len and n_percent <= max_n_percent:
                seen_accessions.add(clean_id)
                rec.id = clean_id
                rec.name = clean_id
                rec.description = ""
                valid_records.append(rec)
                valid_uids.append(uid) # MUST save numeric UID for esummary!
                
            if len(valid_records) >= target_quota:
                break
                
        if len(valid_records) >= target_quota:
            break
            
    if not valid_records:
        print("❌ No sequences passed the QC filters during fetching.", file=sys.stderr)
        sys.exit(1)

    SeqIO.write(valid_records, "results/raw/isolates.fasta", "fasta")
    
    # 4. Fetch Metadata using NUMERIC UIDs (valid_uids)
    print(f"📊 Fetching metadata for the {len(valid_uids)} successful sequences...")
    with Entrez.esummary(db="nucleotide", id=",".join(valid_uids)) as summary_handle:
        summaries = Entrez.read(summary_handle)
    
    meta_data = []
    for s in summaries:
        raw_acc = s.get("Caption") or s.get("AccessionVersion") or ""
        clean_acc = clean_accession(raw_acc)
        
        meta_data.append({
            "Accession": clean_acc,
            "Title": s.get("Title", ""),
            "Length_bp": int(s.get("Length", 0)),
            "UpdateDate": s.get("UpdateDate", "")
        })
            
    pd.DataFrame(meta_data).to_csv("results/raw/metadata.csv", index=False)
    print(f"✅ Success! Saved {len(valid_records)} pristine NCBI sequences and metadata.")

if __name__ == "__main__":
    main()