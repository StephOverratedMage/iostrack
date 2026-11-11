import os
import re
import sys
import pandas as pd

# Allow importing directly from workflow/scripts/
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
    annotate_indel_consequence,
    write_multisample_vcf,
    build_nj_tree,
    tree_to_newick,
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
    ref_seq = "ATGAAATAA"
    cds_models = [
        {
            "gene": "F",
            "contig_id": "test_seg",
            "product": "fusion protein",
            "strand": 1,
            "cds_positions": list(range(9)),
            "pos_to_idx": {i: i for i in range(9)},
            "ref_seq": ref_seq,
        }
    ]

    # 1. Pos 5: AAA (Lys) -> AAG (Lys) => Synonymous
    res_syn = annotate_snp_consequence("test_seg", 5, "A", "G", cds_models)
    assert res_syn["Consequence"] == "Synonymous"
    assert res_syn["AA_Change"] == "F:p.K2="

    # 2. Pos 3: AAA (Lys) -> GAA (Glu) => Missense
    res_mis = annotate_snp_consequence("test_seg", 3, "A", "G", cds_models)
    assert res_mis["Consequence"] == "Missense"
    assert res_mis["AA_Change"] == "F:p.K2E"

    # 3. Pos 3: AAA (Lys) -> TAA (Stop) => Nonsense (Stop)
    res_stop = annotate_snp_consequence("test_seg", 3, "A", "T", cds_models)
    assert res_stop["Consequence"] == "Nonsense (Stop)"
    assert res_stop["AA_Change"] == "F:p.K2*"


def test_annotate_indel_consequence_frameshift_vs_inframe():
    """Verifies that 3-bp indels in CDS are classified as In-Frame and non-3-bp as Frameshift."""
    ref_seq = "ATGAAATAA"
    cds_models = [
        {
            "contig_id": "test_seg",
            "gene": "F",
            "product": "fusion protein",
            "strand": 1,
            "cds_positions": list(range(9)),
            "pos_to_idx": {i: i for i in range(9)},
            "ref_seq": ref_seq,
        }
    ]

    res_fs = annotate_indel_consequence("test_seg", 3, "Deletion", "AA", cds_models)
    assert res_fs["Gene"] == "F"
    assert res_fs["Consequence"] == "Frameshift_Deletion"
    assert res_fs["Indel_Length"] == 2

    res_if = annotate_indel_consequence("test_seg", 3, "Insertion", "GGT", cds_models)
    assert res_if["Gene"] == "F"
    assert res_if["Consequence"] == "In-Frame_Insertion"
    assert res_if["Indel_Length"] == 3

    res_nc = annotate_indel_consequence("test_seg", 50, "Insertion", "A", cds_models)
    assert res_nc["Gene"] == "Noncoding"
    assert res_nc["Consequence"] == "Noncoding_Insertion"


def test_write_multisample_vcf(tmp_path):
    """Verifies 0-based Target_Pos is written as 1-based POS in VCFv4.2 with accurate AF."""
    vcf_out = tmp_path / "test_variants.vcf"
    snp_df = pd.DataFrame(
        [
            {
                "Accession": "ISO1",
                "Target_Contig": "NC_001781",
                "Target_Pos": 99,  # 0-based -> should become 100 in VCF
                "Ref_Base": "A",
                "Alt_Base": "G",
                "Gene": "G",
                "Consequence": "Missense",
                "Codon_Change": "ACA>GCA (pos 1/3)",
                "AA_Change": "G:p.T34A",
            }
        ]
    )
    write_multisample_vcf(str(vcf_out), {"NC_001781": 15225}, snp_df, ["ISO1", "ISO2"])

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


def test_neighbor_joining_tree_construction():
    """Verifies Saitou-Nei Neighbor-Joining groups sister taxa and exports valid Newick."""
    dist_df = pd.DataFrame(
        [
            [0, 2, 7, 8],
            [2, 0, 7, 8],
            [7, 7, 0, 3],
            [8, 8, 3, 0],
        ],
        index=["REF", "ISO_A", "ISO_B", "ISO_C"],
        columns=["REF", "ISO_A", "ISO_B", "ISO_C"],
    )
    tree = build_nj_tree(dist_df, root_label="REF")
    nwk = tree_to_newick(tree)

    assert nwk.endswith(";")
    assert "REF:0.0000" in nwk
    assert "ISO_A" in nwk and "ISO_B" in nwk and "ISO_C" in nwk