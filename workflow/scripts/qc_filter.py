import os
import yaml
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def main():
    # 1. Load thresholds from the config file
    with open("config/config.yaml", "r") as f:
        config = yaml.safe_load(f)
        
    min_len = config["min_length"]
    max_n_pct = config["max_ambiguous_n_percent"]
    
    os.makedirs("results/qc", exist_ok=True)
    
    clean_records = []
    qc_stats = []
    
    # 2. Loop through the raw sequences one by one
    print(f"🔬 Running QC (Min length: {min_len}bp, Max N%: {max_n_pct}%)...")
    for record in SeqIO.parse("results/raw/isolates.fasta", "fasta"):
        seq_str = str(record.seq).upper()
        length = len(seq_str)
        
        # Calculate ambiguous bases (N) and GC content
        n_count = seq_str.count("N")
        n_percent = (n_count / length) * 100 if length > 0 else 100
        gc_percent = gc_fraction(seq_str) * 100
        
        # 3. Apply the filters
        status = "PASS"
        if length < min_len:
            status = "FAIL_LENGTH"
        elif n_percent > max_n_pct:
            status = "FAIL_N_CONTENT"
            
        # Log the stats for our GUI report
        qc_stats.append({
            "Accession": record.id,
            "Length_bp": length,
            "GC_Percent": round(gc_percent, 2),
            "N_Percent": round(n_percent, 2),
            "Status": status
        })
        
        if status == "PASS":
            clean_records.append(record)
            
    # 4. Save the clean sequences and the report
    SeqIO.write(clean_records, "results/qc/clean_isolates.fasta", "fasta")
    
    df = pd.DataFrame(qc_stats)
    df.to_csv("results/qc/qc_report.csv", index=False)
    
    print(f"✅ QC Complete! Kept {len(clean_records)} / {len(qc_stats)} sequences.")
    print("📄 Saved to results/qc/clean_isolates.fasta and qc_report.csv")

if __name__ == "__main__":
    main()