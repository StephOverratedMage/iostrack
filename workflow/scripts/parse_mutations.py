import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
CONSEQUENCE_RANK = {"Nonsense (Stop)": 4, "Missense": 3, "Synonymous": 2, "Non-coding": 1}
INDEL_RANK = {
    "Frameshift_Insertion": 3,
    "Frameshift_Deletion": 3,
    "In-Frame_Insertion": 2,
    "In-Frame_Deletion": 2,
    "Noncoding_Insertion": 1,
    "Noncoding_Deletion": 1,
}


def clean_accession(raw_id):
    """Extracts clean base accession from NCBI, ENA, BV-BRC, or local FASTA headers."""
    cleaned = str(raw_id).strip().strip("|")
    if "|" in cleaned:
        cleaned = cleaned.split("|")[-1]
    match = re.search(r"([A-Z]{1,6}_?[0-9]{4,})", cleaned)
    return match.group(1) if match else cleaned.split(".")[0]


def classify_substitution(ref_base, alt_base):
    """Classifies a SNP as Transition or Transversion."""
    if (ref_base, alt_base) in TRANSITIONS:
        return "Transition"
    return "Transversion"


def load_genbank_cds(gb_path):
    """
    Parses CDS features from reference.gb and builds a fast coordinate-to-codon index.
    Returns (gene_rows_list, cds_models_list) supporting multiple contigs.
    """
    gene_rows = []
    cds_models = []

    if not os.path.exists(gb_path):
        return gene_rows, cds_models

    try:
        gb_records = list(SeqIO.parse(gb_path, "genbank"))
    except Exception:
        return gene_rows, cds_models

    seen_genes = set()

    for gb_record in gb_records:
        contig_id = clean_accession(gb_record.id)
        ref_seq_str = str(gb_record.seq).upper()

        for feat in gb_record.features:
            if feat.type != "CDS":
                continue

            gene_name = (
                feat.qualifiers.get("gene", [None])[0]
                or feat.qualifiers.get("locus_tag", [None])[0]
                or feat.qualifiers.get("product", ["CDS"])[0]
            )
            product_name = feat.qualifiers.get("product", [gene_name])[0]
            codon_start = int(feat.qualifiers.get("codon_start", [1])[0]) - 1

            start_0 = int(feat.location.start)
            end_0 = int(feat.location.end)
            strand_num = feat.location.strand if feat.location.strand is not None else 1
            strand_str = "+" if strand_num >= 0 else "-"

            gene_key = (contig_id, gene_name, start_0 + 1, end_0)
            if gene_key not in seen_genes:
                seen_genes.add(gene_key)
                gene_rows.append({
                    "Contig": contig_id,
                    "Gene": str(gene_name),
                    "Product": str(product_name),
                    "Start": start_0 + 1,
                    "End": end_0,
                    "Strand": strand_str,
                })

            cds_positions = []
            for part in feat.location.parts:
                p_start, p_end = int(part.start), int(part.end)
                if (part.strand or 1) >= 0:
                    cds_positions.extend(range(p_start, p_end))
                else:
                    cds_positions.extend(range(p_end - 1, p_start - 1, -1))

            cds_positions = cds_positions[codon_start:]
            pos_to_idx = {g_pos: idx for idx, g_pos in enumerate(cds_positions)}

            cds_models.append({
                "contig_id": contig_id,
                "gene": str(gene_name),
                "product": str(product_name),
                "strand": strand_num,
                "cds_positions": cds_positions,
                "pos_to_idx": pos_to_idx,
                "ref_seq": ref_seq_str,
            })

    return gene_rows, cds_models


def annotate_snp_consequence(contig_id, pos_0, ref_base, alt_base, cds_models):
    """Determines the Gene, Product, Consequence, Codon_Change, and Amino Acid Change (HGVS p.) for a 0-based SNP."""
    best_hit = {
        "Gene": "Intergenic",
        "Product": "Non-coding region",
        "Consequence": "Non-coding",
        "Codon_Change": "Non-coding",
        "AA_Change": "Non-coding",
    }

    for model in cds_models:
        if model["contig_id"] != contig_id:
            continue

        if pos_0 not in model["pos_to_idx"]:
            continue

        idx = model["pos_to_idx"][pos_0]
        codon_start_idx = (idx // 3) * 3
        codon_g_coords = model["cds_positions"][codon_start_idx : codon_start_idx + 3]

        if len(codon_g_coords) < 3:
            continue

        aa_num = (idx // 3) + 1
        ref_seq = model["ref_seq"]
        strand = model["strand"]

        ref_codon_chars = []
        alt_codon_chars = []
        for g_pos in codon_g_coords:
            r_nt = ref_seq[g_pos] if g_pos < len(ref_seq) else "N"
            a_nt = alt_base if g_pos == pos_0 else r_nt
            if strand < 0:
                r_nt = COMPLEMENT.get(r_nt, "N")
                a_nt = COMPLEMENT.get(a_nt, "N")
            ref_codon_chars.append(r_nt)
            alt_codon_chars.append(a_nt)

        ref_codon_str = "".join(ref_codon_chars)
        alt_codon_str = "".join(alt_codon_chars)
        codon_change_str = f"{ref_codon_str}>{alt_codon_str} (pos {(idx % 3) + 1}/3)"

        ref_aa = str(Seq(ref_codon_str).translate())
        alt_aa = str(Seq(alt_codon_str).translate())

        if ref_aa == alt_aa:
            conseq = "Synonymous"
            aa_str = f"{model['gene']}:p.{ref_aa}{aa_num}="
        elif alt_aa == "*":
            conseq = "Nonsense (Stop)"
            aa_str = f"{model['gene']}:p.{ref_aa}{aa_num}*"
        else:
            conseq = "Missense"
            aa_str = f"{model['gene']}:p.{ref_aa}{aa_num}{alt_aa}"

        if CONSEQUENCE_RANK[conseq] > CONSEQUENCE_RANK[best_hit["Consequence"]]:
            best_hit = {
                "Gene": model["gene"],
                "Product": model["product"],
                "Consequence": conseq,
                "Codon_Change": codon_change_str,
                "AA_Change": aa_str,
            }

    return best_hit


def annotate_indel_consequence(contig_id, pos_0, indel_type, indel_seq, cds_models):
    """
    Classifies an Insertion or Deletion at 0-based pos_0 as Frameshift, In-Frame,
    or Noncoding across overlapping CDS models.
    """
    indel_seq = str(indel_seq).upper()
    indel_len = len(indel_seq)
    short_type = "ins" if indel_type == "Insertion" else "del"

    best_hit = {
        "Gene": "Noncoding",
        "Product": "Non-coding region",
        "Consequence": f"Noncoding_{indel_type}",
        "Indel_Length": indel_len,
        "Protein_Impact": f"g.{pos_0 + 1}_{short_type}{indel_len}bp",
    }

    for model in cds_models:
        if model.get("contig_id") != contig_id:
            continue
        if pos_0 not in model["pos_to_idx"]:
            continue

        cds_idx = model["pos_to_idx"][pos_0]
        aa_pos = (cds_idx // 3) + 1
        gene_name = model["gene"]

        if indel_len % 3 == 0:
            aa_count = indel_len // 3
            conseq = f"In-Frame_{indel_type}"
            impact_str = f"{gene_name}:p.aa{aa_pos}_{short_type}{aa_count}aa ({indel_len}bp)"
        else:
            conseq = f"Frameshift_{indel_type}"
            sign = "+" if indel_type == "Insertion" else "-"
            impact_str = f"{gene_name}:p.aa{aa_pos}fs ({sign}{indel_len}bp)"

        if INDEL_RANK.get(conseq, 0) > INDEL_RANK.get(best_hit["Consequence"], 0):
            best_hit = {
                "Gene": gene_name,
                "Product": model.get("product", ""),
                "Consequence": conseq,
                "Indel_Length": indel_len,
                "Protein_Impact": impact_str,
            }

    return best_hit


def write_multisample_vcf(vcf_path, contig_lengths, snp_df, sample_list):
    """Writes a standards-compliant VCF v4.2 file supporting multiple contigs."""
    with open(vcf_path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=Iostrack_Pathogen_Surveillance\n")

        for c_id, c_len in contig_lengths.items():
            out.write(f"##contig=<ID={c_id},length={c_len}>\n")

        out.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Variant">\n')
        out.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency across mapped cohort">\n')
        out.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Viral Gene Name">\n')
        out.write('##INFO=<ID=CONSEQ,Number=1,Type=String,Description="Functional Consequence">\n')
        out.write('##INFO=<ID=AA,Number=1,Type=String,Description="Amino Acid Change">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (Haploid)">\n')

        sample_cols = "\t".join(sample_list)
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_cols}\n")

        if snp_df.empty or not sample_list:
            return

        total_samples = len(sample_list)
        grouped = snp_df.groupby(["Target_Contig", "Target_Pos", "Ref_Base", "Alt_Base"])

        for (ctg, pos_0, ref_b, alt_b), group in sorted(grouped, key=lambda x: (x[0][0], x[0][1])):
            pos_1 = int(pos_0) + 1
            mut_samples = set(group["Accession"].unique())
            ns = len(mut_samples)
            af = round(ns / total_samples, 4)

            first_row = group.iloc[0]
            gene = str(first_row.get("Gene", "Intergenic")).replace(" ", "_")
            conseq = str(first_row.get("Consequence", "Non-coding")).replace(" ", "_")
            aa = str(first_row.get("AA_Change", "Non-coding")).replace(" ", "_")

            info_str = f"NS={ns};AF={af};GENE={gene};CONSEQ={conseq};AA={aa}"
            gt_list = ["1" if s in mut_samples else "0" for s in sample_list]
            gt_str = "\t".join(gt_list)

            out.write(f"{ctg}\t{pos_1}\t.\t{ref_b}\t{alt_b}\t.\tPASS\t{info_str}\tGT\t{gt_str}\n")

def build_nj_tree(dist_df, root_label=None):
    """
    Builds a rooted phylogenetic tree from a symmetric pairwise distance DataFrame
    using the Saitou & Nei (1987) Neighbor-Joining algorithm.
    """
    labels = list(dist_df.index)
    n = len(labels)
    if n == 0:
        return None
    if n == 1:
        return {"name": labels[0], "length": 0.0, "children": [], "is_leaf": True}

    adj = {i: {} for i in range(n)}
    node_names = {i: labels[i] for i in range(n)}
    active = list(range(n))
    d = {
        (i, j): float(dist_df.iloc[i, j])
        for i in range(n)
        for j in range(n)
    }
    next_node = n
    last_edge = (0, 1, d.get((0, 1), 0.0))

    while len(active) > 2:
        r = len(active)
        row_sums = {i: sum(d[(i, k)] for k in active if k != i) for i in active}

        best_pair = None
        best_q = float("inf")
        for idx_i in range(r):
            i = active[idx_i]
            for idx_j in range(idx_i + 1, r):
                j = active[idx_j]
                q_val = (r - 2) * d[(i, j)] - row_sums[i] - row_sums[j]
                if q_val < best_q:
                    best_q = q_val
                    best_pair = (i, j)

        i, j = best_pair
        dij = d[(i, j)]
        limb_i = 0.5 * dij + (row_sums[i] - row_sums[j]) / (2.0 * (r - 2))
        limb_i = max(0.0, limb_i)
        limb_j = max(0.0, dij - limb_i)

        u = next_node
        next_node += 1
        adj[u] = {}
        node_names[u] = f"Internal_{u}"

        adj[u][i] = limb_i
        adj[i][u] = limb_i
        adj[u][j] = limb_j
        adj[j][u] = limb_j

        for k in active:
            if k not in (i, j):
                duk = max(0.0, 0.5 * (d[(i, k)] + d[(j, k)] - dij))
                d[(u, k)] = duk
                d[(k, u)] = duk
                d[(u, u)] = 0.0

        active.remove(i)
        active.remove(j)
        active.append(u)

    if len(active) == 2:
        u1, u2 = active[0], active[1]
        rem_dist = max(0.0, d.get((u1, u2), 0.0))
        adj[u1][u2] = rem_dist
        adj[u2][u1] = rem_dist
        last_edge = (u1, u2, rem_dist)

    def build_directed(curr, parent, branch_len):
        children = []
        for nbr, w in adj[curr].items():
            if nbr != parent:
                children.append(build_directed(nbr, curr, w))
        children.sort(key=lambda c: (0 if c["is_leaf"] else 1, c["sub_size"], c["name"]))
        sub_size = 1 if not children else sum(c["sub_size"] for c in children)
        return {
            "name": node_names.get(curr, str(curr)),
            "length": float(branch_len),
            "children": children,
            "is_leaf": len(children) == 0,
            "sub_size": sub_size,
        }

    if root_label and root_label in labels:
        root_idx = labels.index(root_label)
        nbr_idx, w = next(iter(adj[root_idx].items()))
        outgroup_node = {
            "name": root_label,
            "length": 0.0,
            "children": [],
            "is_leaf": True,
            "sub_size": 1,
        }
        ingroup_node = build_directed(nbr_idx, root_idx, w)
        return {
            "name": "Root",
            "length": 0.0,
            "children": [outgroup_node, ingroup_node],
            "is_leaf": False,
            "sub_size": n,
        }

    u1, u2, rem_dist = last_edge
    left_sub = build_directed(u1, u2, rem_dist / 2.0)
    right_sub = build_directed(u2, u1, rem_dist / 2.0)
    return {
        "name": "Root",
        "length": 0.0,
        "children": [left_sub, right_sub],
        "is_leaf": False,
        "sub_size": n,
    }


def tree_to_newick(tree_node):
    """Serializes a rooted tree dictionary into standard Newick (.nwk) format."""
    if tree_node is None:
        return ";"

    def _recurse(node):
        if not node["children"]:
            safe = re.sub(r"[:(),;\s]+", "_", str(node["name"]))
            return f"{safe}:{node['length']:.4f}"
        inner = ",".join(_recurse(c) for c in node["children"])
        return f"({inner}):{node['length']:.4f}"

    if not tree_node["children"]:
        safe = re.sub(r"[:(),;\s]+", "_", str(tree_node["name"]))
        return f"({safe}:{tree_node['length']:.4f});"
    return "(" + ",".join(_recurse(c) for c in tree_node["children"]) + ");"

def main():
    os.makedirs("results/summary", exist_ok=True)

    ref_records = list(SeqIO.parse("results/raw/reference.fasta", "fasta"))
    contig_lengths = {clean_accession(rec.id): len(rec.seq) for rec in ref_records}

    gene_rows, cds_models = load_genbank_cds("results/raw/reference.gb")
    genes_df = pd.DataFrame(gene_rows) if gene_rows else pd.DataFrame(columns=["Contig", "Gene", "Product", "Start", "End", "Strand"])
    genes_df.to_csv("results/summary/gene_features.csv", index=False)

    meta_df = pd.read_csv("results/raw/metadata.csv")
    qc_df = pd.read_csv("results/qc/qc_report.csv")

    meta_df["Accession"] = meta_df["Accession"].astype(str).apply(clean_accession)
    qc_df["Accession"] = qc_df["Accession"].astype(str).apply(clean_accession)
    meta_df = meta_df.drop_duplicates(subset=["Accession"]).drop(columns=["Length_bp"], errors="ignore")
    qc_df = qc_df.drop_duplicates(subset=["Accession"])

    merged = pd.merge(qc_df, meta_df, on="Accession", how="left")

    mutation_details = []
    with open("results/align/alignments.paf", "r") as f:
        for line in f:
            if "tp:A:P" not in line and "tp:A:I" not in line:
                continue

            parts = line.strip().split("\t")
            isolate_id = clean_accession(parts[0])
            target_contig = clean_accession(parts[5])

            cs_tag = ""
            for tag in parts[12:]:
                if tag.startswith("cs:Z:"):
                    cs_tag = tag.replace("cs:Z:", "")
                    break

            mutation_details.append({
                "Accession": isolate_id,
                "Target_Contig": target_contig,
                "Target_Start": int(parts[7]),
                "Target_End": int(parts[8]),
                "Query_Start": int(parts[2]),
                "Query_End": int(parts[3]),
                "Strand": parts[4],
                "cs_tag": cs_tag,
            })

    details_df = pd.DataFrame(mutation_details)

    indel_cols = [
        "Accession", "Target_Contig", "Target_Pos", "Position_1Based",
        "Indel_Type", "Length_bp", "Sequence", "Gene", "Product",
        "Consequence", "Protein_Impact",
    ]

    if not details_df.empty:
        details_df["Accession"] = details_df["Accession"].astype(str).str.strip()
        details_df["Align_Length"] = details_df["Target_End"] - details_df["Target_Start"]

        snp_coordinates = []
        indel_coordinates = []
        final_mutations = []

        for _, row in details_df.iterrows():
            current_pos = row["Target_Start"]
            target_ctg = row["Target_Contig"]

            snp_count, insertion_count, deletion_count = 0, 0, 0
            ti_count, tv_count = 0, 0
            missense_count, synonymous_count, noncoding_count = 0, 0, 0
            frameshift_count, inframe_count = 0, 0

            tokens = re.findall(r"(:[0-9]+|\*[a-zA-Z]{2}|\+[a-zA-Z]+|-[a-zA-Z]+)", row["cs_tag"])

            for token in tokens:
                op = token[0]
                if op == ":":
                    current_pos += int(token[1:])
                elif op == "*":
                    ref_b = token[1].upper()
                    alt_b = token[2].upper()

                    if ref_b in "ACGT" and alt_b in "ACGT":
                        mut_class = classify_substitution(ref_b, alt_b)
                        if mut_class == "Transition":
                            ti_count += 1
                        else:
                            tv_count += 1

                        ann = annotate_snp_consequence(target_ctg, current_pos, ref_b, alt_b, cds_models)
                        if ann["Consequence"] in ("Missense", "Nonsense (Stop)"):
                            missense_count += 1
                        elif ann["Consequence"] == "Synonymous":
                            synonymous_count += 1
                        else:
                            noncoding_count += 1

                        snp_coordinates.append({
                            "Accession": row["Accession"],
                            "Target_Contig": target_ctg,
                            "Target_Pos": current_pos,
                            "Ref_Base": ref_b,
                            "Alt_Base": alt_b,
                            "Mutation": f"{ref_b}>{alt_b}",
                            "Mut_Class": mut_class,
                            "Gene": ann["Gene"],
                            "Product": ann["Product"],
                            "Consequence": ann["Consequence"],
                            "Codon_Change": ann["Codon_Change"],
                            "AA_Change": ann["AA_Change"],
                        })
                        snp_count += 1
                    current_pos += 1
                elif op == "-":
                    del_seq = token[1:].upper()
                    indel_ann = annotate_indel_consequence(target_ctg, current_pos, "Deletion", del_seq, cds_models)
                    if "Frameshift" in indel_ann["Consequence"]:
                        frameshift_count += 1
                    elif "In-Frame" in indel_ann["Consequence"]:
                        inframe_count += 1

                    indel_coordinates.append({
                        "Accession": row["Accession"],
                        "Target_Contig": target_ctg,
                        "Target_Pos": current_pos,
                        "Position_1Based": current_pos + 1,
                        "Indel_Type": "Deletion",
                        "Length_bp": len(del_seq),
                        "Sequence": del_seq if len(del_seq) <= 30 else del_seq[:27] + "...",
                        "Gene": indel_ann["Gene"],
                        "Product": indel_ann["Product"],
                        "Consequence": indel_ann["Consequence"],
                        "Protein_Impact": indel_ann["Protein_Impact"],
                    })
                    current_pos += len(del_seq)
                    deletion_count += 1
                elif op == "+":
                    ins_seq = token[1:].upper()
                    indel_ann = annotate_indel_consequence(target_ctg, current_pos, "Insertion", ins_seq, cds_models)
                    if "Frameshift" in indel_ann["Consequence"]:
                        frameshift_count += 1
                    elif "In-Frame" in indel_ann["Consequence"]:
                        inframe_count += 1

                    indel_coordinates.append({
                        "Accession": row["Accession"],
                        "Target_Contig": target_ctg,
                        "Target_Pos": current_pos,
                        "Position_1Based": current_pos + 1,
                        "Indel_Type": "Insertion",
                        "Length_bp": len(ins_seq),
                        "Sequence": ins_seq if len(ins_seq) <= 30 else ins_seq[:27] + "...",
                        "Gene": indel_ann["Gene"],
                        "Product": indel_ann["Product"],
                        "Consequence": indel_ann["Consequence"],
                        "Protein_Impact": indel_ann["Protein_Impact"],
                    })
                    insertion_count += 1

            segment_titv = round(ti_count / tv_count, 2) if tv_count > 0 else float(ti_count)

            final_mutations.append({
                "Accession": row["Accession"],
                "Target_Contig": target_ctg,
                "Target_Start": row["Target_Start"],
                "Target_End": row["Target_End"],
                "Query_Start": row["Query_Start"],
                "Query_End": row["Query_End"],
                "Strand": row["Strand"],
                "SNPs": snp_count,
                "Transitions": ti_count,
                "Transversions": tv_count,
                "Ti_Tv_Ratio": segment_titv,
                "Missense_SNPs": missense_count,
                "Synonymous_SNPs": synonymous_count,
                "Noncoding_SNPs": noncoding_count,
                "Insertions": insertion_count,
                "Deletions": deletion_count,
                "Frameshift_Indels": frameshift_count,
                "InFrame_Indels": inframe_count,
            })

        clean_details_df = pd.DataFrame(final_mutations)

        detailed_full = pd.merge(merged, clean_details_df, on="Accession", how="inner")
        detailed_full.to_csv("results/summary/split_alignments.csv", index=False)

        sum_cols = [
            "SNPs", "Transitions", "Transversions",
            "Missense_SNPs", "Synonymous_SNPs", "Noncoding_SNPs",
            "Insertions", "Deletions", "Frameshift_Indels", "InFrame_Indels",
        ]
        mut_agg = clean_details_df.groupby("Accession", as_index=False)[sum_cols].sum()

        mut_agg["Ti_Tv_Ratio"] = mut_agg.apply(
            lambda r: round(r["Transitions"] / r["Transversions"], 2) if r["Transversions"] > 0 else float(r["Transitions"]),
            axis=1,
        )
        mut_agg["Missense_Syn_Ratio"] = mut_agg.apply(
            lambda r: round(r["Missense_SNPs"] / r["Synonymous_SNPs"], 2) if r["Synonymous_SNPs"] > 0 else float(r["Missense_SNPs"]),
            axis=1,
        )

        final_df = pd.merge(merged, mut_agg, on="Accession", how="inner")

        snp_df = pd.DataFrame(snp_coordinates)
        if snp_df.empty:
            snp_df = pd.DataFrame(columns=[
                "Accession", "Target_Contig", "Target_Pos", "Ref_Base", "Alt_Base", "Mutation",
                "Mut_Class", "Gene", "Product", "Consequence", "Codon_Change", "AA_Change",
            ])
        snp_df.to_csv("results/summary/snp_coordinates.csv", index=False)

        indel_df = pd.DataFrame(indel_coordinates, columns=indel_cols) if indel_coordinates else pd.DataFrame(columns=indel_cols)
        indel_df.to_csv("results/summary/indel_coordinates.csv", index=False)
    else:
        final_df = merged.iloc[0:0].copy()
        for col in [
            "SNPs", "Transitions", "Transversions", "Ti_Tv_Ratio",
            "Missense_SNPs", "Synonymous_SNPs", "Noncoding_SNPs", "Missense_Syn_Ratio",
            "Insertions", "Deletions", "Frameshift_Indels", "InFrame_Indels",
        ]:
            final_df[col] = 0
        pd.DataFrame(columns=[
            "Accession", "Target_Contig", "Target_Start", "Target_End", "Query_Start", "Query_End",
            "Strand", "SNPs", "Transitions", "Transversions", "Ti_Tv_Ratio",
            "Missense_SNPs", "Synonymous_SNPs", "Noncoding_SNPs", "Insertions", "Deletions",
            "Frameshift_Indels", "InFrame_Indels",
        ]).to_csv("results/summary/split_alignments.csv", index=False)
        snp_df = pd.DataFrame(columns=[
            "Accession", "Target_Contig", "Target_Pos", "Ref_Base", "Alt_Base", "Mutation",
            "Mut_Class", "Gene", "Product", "Consequence", "Codon_Change", "AA_Change",
        ])
        snp_df.to_csv("results/summary/snp_coordinates.csv", index=False)
        pd.DataFrame(columns=indel_cols).to_csv("results/summary/indel_coordinates.csv", index=False)

    final_df = final_df[final_df["Status"] == "PASS"]
    int_cols = [
        "SNPs", "Transitions", "Transversions",
        "Missense_SNPs", "Synonymous_SNPs", "Noncoding_SNPs",
        "Insertions", "Deletions", "Frameshift_Indels", "InFrame_Indels",
    ]
    final_df[int_cols] = final_df[int_cols].fillna(0).astype(int)
    final_df["Ti_Tv_Ratio"] = final_df["Ti_Tv_Ratio"].fillna(0.0)
    final_df["Missense_Syn_Ratio"] = final_df["Missense_Syn_Ratio"].fillna(0.0)
    final_df.to_csv("results/summary/sample_report.csv", index=False)

    mapped_samples = final_df["Accession"].tolist()
    write_multisample_vcf("results/summary/variants.vcf", contig_lengths, snp_df, mapped_samples)


if __name__ == "__main__":
    main()