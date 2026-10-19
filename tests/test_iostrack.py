import os
import re
import sys
import pandas as pd

# Allow importing directly from scripts/
sys.path.insert(
    0,
    os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts")
    ),
)

from fetch_ena import clean_accession
from parse_mutations import (
    classify_substitution,
    annotate_snp_consequence,
    write_multisample_vcf,
)


def test_clean_accession_formats():
    """Verifies NCBI versioned accessions and ENA piped headers strip cleanly."""
    assert clean_accession("NC_001781.1") == "NC_001781"
    assert clean_accession("ENA|OR287049|OR287049.1") == "OR287049"
    assert clean_accession("gi|9629198|ref|NC_001781.1|") == "NC_001781"


def test_transition_transversion_classification():
    """Verifies canonical purine/pyrimidine transitions vs. transversions."""
    for ref_b, alt_b in [("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")]:
        assert classify_substitution(ref_b, alt_b) == "Transition"

    for ref_b, alt_b in [("A", "T"), ("A", "C"), ("G", "T"), ("C", "A")]:
        assert classify_substitution(ref_b, alt_b) == "Transversion"


def test_cs_tag_coordinate_state_machine():
    """
    Ensures deletions (-acg) advance reference coordinates by len(del),
    while insertions (+tttt) do NOT advance reference coordinates.
    """
    # Start at 0-based pos 100:
    # :10 -> pos 110
    # *ag -> SNP at 110, pos becomes 111
    # +tttt -> insertion in query, ref pos MUST stay 111
    # -acg -> 3bp deletion in query, ref pos advances by 3 -> 114
    # *ct -> SNP at 114, pos becomes 115
    cs_tag = ":10*ag+tttt-acg*ct"
    current_pos = 100
    called_snps = []

    tokens = re.findall(r"(:[0-9]+|\*[a-zA-Z]{2}|\+[a-zA-Z]+|-[a-zA-Z]+)", cs_tag)
    for token in tokens:
        op = token[0]
        if op == ":":
            current_pos += int(token[1:])
        elif op == "*":
            called_snps.append((current_pos, token[1].upper(), token[2].upper()))
            current_pos += 1
        elif op == "-":
            current_pos += len(token) - 1
        elif op == "+":
            pass  # Insertions do not consume reference bases

    assert called_snps == [(110, "A", "G"), (114, "C", "T")]
    assert current_pos == 115


def test_annotate_snp_consequence_coding_and_stop():
    """Tests codon translation for Synonymous, Missense, and Nonsense SNPs."""
    # Toy reference: ATG (Met, aa 1) + AAA (Lys, aa 2) + TAA (Stop, aa 3)
    ref_seq = "ATGAAATAA"
    cds_models = [
        {
            "gene": "F",
            "product": "fusion protein",
            "strand": 1,
            "cds_positions": list(range(9)),
            "pos_to_idx": {i: i for i in range(9)},
            "ref_seq": ref_seq,
        }
    ]

    # 1. Pos 5: AAA (Lys) -> AAG (Lys) => Synonymous
    res_syn = annotate_snp_consequence(5, "A", "G", cds_models)
    assert res_syn["Consequence"] == "Synonymous"
    assert res_syn["AA_Change"] == "F:p.K2="

    # 2. Pos 3: AAA (Lys) -> GAA (Glu) => Missense
    res_mis = annotate_snp_consequence(3, "A", "G", cds_models)
    assert res_mis["Consequence"] == "Missense"
    assert res_mis["AA_Change"] == "F:p.K2E"

    # 3. Pos 3: AAA (Lys) -> TAA (Stop) => Nonsense (Stop)
    res_stop = annotate_snp_consequence(3, "A", "T", cds_models)
    assert res_stop["Consequence"] == "Nonsense (Stop)"
    assert res_stop["AA_Change"] == "F:p.K2*"


def test_write_multisample_vcf(tmp_path):
    """Verifies 0-based Target_Pos is written as 1-based POS in VCFv4.2 with accurate AF."""
    vcf_out = tmp_path / "test_variants.vcf"
    snp_df = pd.DataFrame(
        [
            {
                "Accession": "ISO1",
                "Target_Pos": 99,  # 0-based -> should become 100 in VCF
                "Ref_Base": "A",
                "Alt_Base": "G",
                "Gene": "G",
                "Consequence": "Missense",
                "AA_Change": "G:p.T34A",
            }
        ]
    )
    write_multisample_vcf(str(vcf_out), "NC_001781", 15225, snp_df, ["ISO1", "ISO2"])

    lines = vcf_out.read_text().strip().splitlines()
    data_lines = [l for l in lines if not l.startswith("#")]
    assert len(data_lines) == 1

    fields = data_lines[0].split("\t")
    assert fields[0] == "NC_001781"
    assert fields[1] == "100"  # 1-based POS
    assert fields[3] == "A"
    assert fields[4] == "G"
    assert "NS=1;AF=0.5;GENE=G;CONSEQ=Missense;AA=G:p.T34A" in fields[7]
    assert fields[9:] == ["1", "0"]