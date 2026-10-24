import os
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# Import your pipeline's actual functions
sys.path.insert(0, "workflow/scripts" if os.path.exists("workflow/scripts") else "scripts")
import parse_mutations as pm


def run_verification():
    print("🔬 STARTING IOSTRACK MATHEMATICAL & SCIENTIFIC AUDIT...\n")
    errors = 0

    # --- TEST 1: SUBSTITUTION & CODON TRANSLATION LOGIC ---
    print("1️⃣ Testing Biochemical Substitution & Codon Translation Logic...")
    assert pm.classify_substitution("A", "G") == "Transition", "Failed A>G Transition"
    assert pm.classify_substitution("C", "T") == "Transition", "Failed C>T Transition"
    assert pm.classify_substitution("A", "T") == "Transversion", "Failed A>T Transversion"

    # Build a synthetic 9-bp gene: ATG (Met/M) - GAA (Glu/E) - TTT (Phe/F)
    synthetic_model = [{
        "contig_id": "TEST_CTG",
        "gene": "TEST_GENE",
        "product": "Test Protein",
        "strand": 1,
        "cds_positions": [0, 1, 2, 3, 4, 5, 6, 7, 8],
        "pos_to_idx": {i: i for i in range(9)},
        "ref_seq": "ATGGAATTT"
    }]

    # 1a. Synonymous: GAA (Glu) -> GAG (Glu) at 0-based pos 5
    res_syn = pm.annotate_snp_consequence("TEST_CTG", 5, "A", "G", synthetic_model)
    assert res_syn["Consequence"] == "Synonymous" and res_syn["AA_Change"] == "TEST_GENE:p.E2=", f"Failed Synonymous: {res_syn}"

    # 1b. Missense: GAA (Glu) -> GTA (Val) at 0-based pos 4
    res_mis = pm.annotate_snp_consequence("TEST_CTG", 4, "A", "T", synthetic_model)
    assert res_mis["Consequence"] == "Missense" and res_mis["AA_Change"] == "TEST_GENE:p.E2V", f"Failed Missense: {res_mis}"

    # 1c. Nonsense (Stop): GAA (Glu) -> TAA (Stop/*) at 0-based pos 3
    res_stop = pm.annotate_snp_consequence("TEST_CTG", 3, "G", "T", synthetic_model)
    assert res_stop["Consequence"] == "Nonsense (Stop)" and res_stop["AA_Change"] == "TEST_GENE:p.E2*", f"Failed Nonsense: {res_stop}"
    print("   ✅ Passed: Synonymous (p.E2=), Missense (p.E2V), and Stop-gain (p.E2*) math is exact.\n")

    # --- TEST 2: AUDIT CURRENT PIPELINE OUTPUTS ---
    report_path = "results/summary/sample_report.csv"
    snps_path = "results/summary/snp_coordinates.csv"
    genes_path = "results/summary/gene_features.csv"

    if not os.path.exists(report_path):
        print("ℹ️ No active 'results/summary/sample_report.csv' found to audit. Run a pipeline query first.")
        return

    print("2️⃣ Auditing Current Pipeline Output Tables for Algebraic Conservation...")
    df = pd.read_csv(report_path)
    snps_df = pd.read_csv(snps_path) if os.path.exists(snps_path) else pd.DataFrame()
    genes_df = pd.read_csv(genes_path) if os.path.exists(genes_path) else pd.DataFrame()

    for _, row in df.iterrows():
        acc = row["Accession"]
        snps = int(row["SNPs"])
        ti, tv = int(row["Transitions"]), int(row["Transversions"])
        mis, syn, non = int(row["Missense_SNPs"]), int(row["Synonymous_SNPs"]), int(row["Noncoding_SNPs"])

        if snps != (ti + tv):
            print(f"   ❌ ERROR in {acc}: SNPs ({snps}) != Ti ({ti}) + Tv ({tv})")
            errors += 1
        if snps != (mis + syn + non):
            print(f"   ❌ ERROR in {acc}: SNPs ({snps}) != Missense ({mis}) + Syn ({syn}) + Noncoding ({non})")
            errors += 1

        # Cross-check against raw snp_coordinates.csv rows
        if not snps_df.empty:
            coord_count = len(snps_df[snps_df["Accession"] == acc])
            if snps != coord_count:
                print(f"   ❌ ERROR in {acc}: sample_report SNPs ({snps}) != snp_coordinates rows ({coord_count})")
                errors += 1

    if errors == 0:
        print(f"   ✅ Passed: All {len(df)} mapped isolates obey strict SNP, Ti/Tv, and CDS partition equations.\n")

    # --- TEST 3: CONTIG & COORDINATE BOUNDARY AUDIT ---
    print("3️⃣ Auditing Reference Contig & GenBank CDS Coordinate Boundaries...")
    if not snps_df.empty and not genes_df.empty and "Target_Contig" in snps_df.columns:
        ref_records = {pm.clean_accession(r.id): len(r.seq) for r in SeqIO.parse("results/raw/reference.fasta", "fasta")}
        for _, s_row in snps_df.iterrows():
            ctg = s_row["Target_Contig"]
            pos = int(s_row["Target_Pos"])
            if ctg not in ref_records:
                print(f"   ❌ ERROR: SNP mapped to unknown contig {ctg}")
                errors += 1
            elif not (0 <= pos < ref_records[ctg]):
                print(f"   ❌ ERROR: SNP position {pos} out of bounds for {ctg} (length {ref_records[ctg]})")
                errors += 1
        if errors == 0:
            print(f"   ✅ Passed: All {len(snps_df)} SNP coordinates fall strictly within their Reference Contig bounds.\n")

    if errors == 0:
        print("🏆 AUDIT COMPLETE: 0 mathematical or coordinate discrepancies detected!")
    else:
        print(f"⚠️ AUDIT FAILED with {errors} error(s).")


if __name__ == "__main__":
    run_verification()