import os
import shutil
import subprocess
import pandas as pd
import plotly.express as px
import streamlit as st
import yaml
import re


# --- THEME SWITCHER LOGIC ---
def set_theme(mode):
    os.makedirs(".streamlit", exist_ok=True)
    if mode == "Dark":
        theme_str = '[theme]\nbase="dark"\nprimaryColor="#00B4D8"\nbackgroundColor="#0E1117"\nsecondaryBackgroundColor="#262730"\ntextColor="#FAFAFA"\nfont="sans serif"'
    else:
        theme_str = '[theme]\nbase="light"\nprimaryColor="#00B4D8"\nbackgroundColor="#FFFFFF"\nsecondaryBackgroundColor="#F0F2F6"\ntextColor="#31333F"\nfont="sans serif"'
    
    with open(".streamlit/config.toml", "w") as f:
        f.write(theme_str)
    st.rerun()

# Determine the current theme by reading the file
current_mode = "Dark"
if os.path.exists(".streamlit/config.toml"):
    with open(".streamlit/config.toml", "r") as f:
        if 'base="light"' in f.read().replace(' ', ''):
            current_mode = "Light"
# ----------------------------

# Set up the webpage layout
st.set_page_config(page_title="Iostrack Dashboard", layout="wide")
st.title("🧬 Iostrack: Pathogen Surveillance Dashboard")

# Place the theme buttons at the very top of the sidebar
st.sidebar.markdown("### 🌓 Appearance")
theme_col1, theme_col2 = st.sidebar.columns(2)
if theme_col1.button("🌙 Dark", width="stretch"):
    if current_mode != "Dark": 
        set_theme("Dark")
if theme_col2.button("☀️ Light", width="stretch"):
    if current_mode != "Light": 
        set_theme("Light")
st.sidebar.markdown("---")

# 1. Read current config to populate the sidebar
with open("config/config.yaml", "r") as f:
    config = yaml.safe_load(f)

# 2. Build the interactive sidebar for config tuning
st.sidebar.header("⚙️ Pipeline Configuration")

data_source = st.sidebar.selectbox(
    "Data Source", 
    ["NCBI", "ENA", "BV-BRC"], 
    index=["NCBI", "ENA", "BV-BRC"].index(config.get("data_source", "NCBI"))
)
config["data_source"] = data_source

@st.cache_data
def load_pathogen_db(filepath="config/pathogens.csv"):
    curated = pd.DataFrame([
        {"Organism": "Severe acute respiratory syndrome coronavirus 2", "Taxon_ID": 2697049, "Default_Ref": "NC_045512", "Default_Min_Len": 28000, "Aliases": "SARS-CoV-2 COVID-19 Betacoronavirus pandemicum"},
        {"Organism": "Respiratory syncytial virus", "Taxon_ID": 11250, "Default_Ref": "NC_001781", "Default_Min_Len": 14000, "Aliases": "RSV Orthopneumovirus hominis Human orthopneumovirus"},
        {"Organism": "Alphainfluenzavirus influenzae", "Taxon_ID": 11320, "Default_Ref": "GCF_000865725.1", "Default_Min_Len": 800, "Aliases": "Influenza A virus Flu A H1N1 H3N2 H5N1"},
        {"Organism": "Betainfluenzavirus influenzae", "Taxon_ID": 11520, "Default_Ref": "GCF_000820495.2", "Default_Min_Len": 800, "Aliases": "Influenza B virus Flu B"},
        {"Organism": "Orthopoxvirus monkeypox", "Taxon_ID": 10244, "Default_Ref": "NC_063383", "Default_Min_Len": 45000, "Aliases": "Monkeypox virus Mpox MPXV"},
        {"Organism": "Orthoflavivirus zikaense", "Taxon_ID": 64320, "Default_Ref": "NC_012532", "Default_Min_Len": 10000, "Aliases": "Zika virus ZIKV"},
        {"Organism": "Orthoflavivirus denguei", "Taxon_ID": 12637, "Default_Ref": "NC_001477", "Default_Min_Len": 10000, "Aliases": "Dengue virus DENV"},
        {"Organism": "Orthoebolavirus zairense", "Taxon_ID": 186538, "Default_Ref": "NC_002549", "Default_Min_Len": 18000, "Aliases": "Zaire ebolavirus Ebola virus EBOV"},
        {"Organism": "Morbillivirus hominis", "Taxon_ID": 11234, "Default_Ref": "NC_001498", "Default_Min_Len": 15000, "Aliases": "Measles morbillivirus Measles virus"},
        {"Organism": "Human immunodeficiency virus 1", "Taxon_ID": 11676, "Default_Ref": "NC_001802", "Default_Min_Len": 8500, "Aliases": "HIV-1 Lentivirus humimdef1"},
    ])
    if os.path.exists(filepath):
        ext_df = pd.read_csv(filepath)
        for col in ["Default_Ref", "Default_Min_Len", "Aliases"]:
            if col not in ext_df.columns:
                ext_df[col] = ""
        combined = pd.concat([curated, ext_df], ignore_index=True)
        return combined.drop_duplicates(subset=["Organism"], keep="first")
    return curated

pathogen_df = load_pathogen_db()

# Initialize defaults so no variable is ever undefined
new_term = config.get("search_term", "")
ena_query = config.get("ena_query", "tax_tree(11250)")
bvbrc_query = config.get("bvbrc_query", "eq(taxon_id,11250)")
new_email = config.get("ncbi_email", "")
auto_ref = config.get("reference_accession", "NC_001781")
auto_tax = str(config.get("ena_taxon", "11250"))
auto_len = int(config.get("min_length", 10000))

# Universal Search Mode (works for both NCBI and ENA)
saved_mode = config.get("search_mode", "By Organism & Keyword")
search_mode = st.sidebar.radio(
    "Search Mode",
    ["By Organism & Keyword", "By Accession List"],
    index=0 if saved_mode == "By Organism & Keyword" else 1
)
config["search_mode"] = search_mode

if search_mode == "By Organism & Keyword":
    st.sidebar.markdown("### 1. Organism Selection")
    search_query = st.sidebar.text_input(
        "🔍 Search Database", 
        value="", 
        placeholder="Type name or synonym (e.g., 'SARS', 'severe acute', 'RSV')"
    )

    if search_query.strip():
        q = search_query.strip()
        mask = (
            pathogen_df["Organism"].astype(str).str.contains(q, case=False, na=False) |
            pathogen_df["Aliases"].astype(str).str.contains(q, case=False, na=False) |
            pathogen_df["Taxon_ID"].astype(str).str.contains(q, case=False, na=False)
        )
        filtered_df = pathogen_df[mask]
        organism_options = filtered_df["Organism"].tolist()[:100]
        
        if not organism_options:
            st.sidebar.warning(f"⚠️ No preset found matching **'{q}'**. Using it as a custom organism query.")
            organism_options = [q]
    else:
        organism_options = pathogen_df["Organism"].head(10).tolist()
        saved_org = config.get("organism", "Respiratory syncytial virus")
        if saved_org and saved_org not in organism_options:
            organism_options = [saved_org] + organism_options

    org_val = st.sidebar.selectbox("Select Organism", options=organism_options) or search_query.strip() or "Respiratory syncytial virus"

    selected_row = pathogen_df[pathogen_df["Organism"] == org_val]
    if not selected_row.empty:
        raw_ref = selected_row.iloc[0].get("Default_Ref", None)
        auto_ref = str(raw_ref) if pd.notna(raw_ref) and str(raw_ref).strip() else auto_ref
        raw_tax = selected_row.iloc[0].get("Taxon_ID", "11250")
        auto_tax = str(int(float(raw_tax))) if pd.notna(raw_tax) and str(raw_tax).strip() else "11250"
        raw_len = selected_row.iloc[0].get("Default_Min_Len", None)
        auto_len = int(float(raw_len)) if pd.notna(raw_len) and str(raw_len).strip() else auto_len

    if st.session_state.get("last_selected_org") != org_val:
        st.session_state["last_selected_org"] = org_val
        if not selected_row.empty and pd.notna(selected_row.iloc[0].get("Default_Ref")) and str(selected_row.iloc[0].get("Default_Ref")).strip():
            st.session_state["ref_acc_input"] = auto_ref

    config["organism"] = org_val
    config["ena_taxon"] = auto_tax

    if data_source == "NCBI":
        type_val = st.sidebar.text_input("Sequence Type (Title Filter)", config.get("genome_type", ""))
        clean_org = str(org_val or "").strip().strip('"')
        clean_type = str(type_val or "").strip().strip('"')
        
        if clean_type:
            new_term = f'"{clean_org}"[Organism] AND "{clean_type}"[Title]'
        else:
            new_term = f'"{clean_org}"[Organism]'
            
        config["genome_type"] = type_val
        config["search_term"] = new_term
        st.sidebar.info(f"**NCBI Query:**\n`{new_term}`")
        
    elif data_source == "ENA":
        taxon_id = st.sidebar.text_input("Taxon ID (auto-filled from Organism)", value=auto_tax)
        ena_query = f"tax_tree({taxon_id.strip()})"
        config["ena_taxon"] = taxon_id.strip()
        config["ena_query"] = ena_query
        st.sidebar.info(f"**ENA Query:**\n`{ena_query}`")
        
    elif data_source == "BV-BRC":
        taxon_id = st.sidebar.text_input("Taxon ID (auto-filled from Organism)", value=auto_tax)
        bvbrc_query = f"eq(taxon_id,{taxon_id.strip()})"
        config["ena_taxon"] = taxon_id.strip() # We can reuse the taxon parameter
        config["bvbrc_query"] = bvbrc_query
        st.sidebar.info(f"**BV-BRC Query:**\n`{bvbrc_query}`")

else:
    # --- BY ACCESSION LIST MODE (NCBI & ENA) ---
    st.sidebar.markdown("### 1. Wild Isolate Accessions")
    acc_val = st.sidebar.text_area(
        "Accessions (comma, space, or newline separated)",
        config.get("accession_list", "OR287049, OR287050"),
        help="Enter NCBI/ENA nucleotide accessions (e.g. OR287049, MT291826) or assembly accessions (e.g. GCA_020736045.1)."
    )
    config["accession_list"] = acc_val
    raw_accs = [x.strip() for x in re.split(r"[,\s]+", str(acc_val)) if x.strip()]

    if data_source == "NCBI":
        new_term = " OR ".join([f"{a}[Accession]" for a in raw_accs]) if raw_accs else ""
        config["search_term"] = new_term
        st.sidebar.info(f"**NCBI Accession Query ({len(raw_accs)}):**\n`{new_term}`")
    elif data_source == "ENA":
        ena_accs = [re.sub(r"\.[0-9]+$", "", a) for a in raw_accs]
        ena_query = "(" + " OR ".join([f'accession="{a}"' for a in ena_accs]) + ")" if ena_accs else ""
        config["ena_query"] = ena_query
        st.sidebar.info(f"**ENA Accession Query ({len(ena_accs)}):**\n`{ena_query}`")
        
    elif data_source == "BV-BRC":
        bvbrc_accs = [re.sub(r"\.[0-9]+$", "", a) for a in raw_accs]
        bvbrc_query = "in(accession,(" + ",".join(bvbrc_accs) + "))" if bvbrc_accs else ""
        config["bvbrc_query"] = bvbrc_query
        st.sidebar.info(f"**BV-BRC Accession Query ({len(bvbrc_accs)}):**\n`{bvbrc_query}`")

# --- BASELINE REFERENCE SECTION ---
st.sidebar.markdown("### 2. Baseline Reference")
if "ref_acc_input" not in st.session_state:
    st.session_state["ref_acc_input"] = config.get("reference_accession", auto_ref)

if search_mode == "By Organism & Keyword":
    if st.sidebar.button("🔍 Auto-Find Reference for Organism", width="stretch"):
        from Bio import Entrez
        Entrez.email = config.get("ncbi_email") or "test@example.com"
        try:
            with Entrez.esearch(db="assembly", term=f'"{config["organism"]}"[Organism] AND "latest refseq"[Filter]', retmax=5) as h:
                asm_res = Entrez.read(h)
            if asm_res["IdList"]:
                with Entrez.esummary(db="assembly", id=asm_res["IdList"][0]) as sh:
                    asm_sum = Entrez.read(sh)
                found_acc = asm_sum["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]
                st.session_state["ref_acc_input"] = found_acc
                st.sidebar.success(f"Found RefSeq Assembly: {found_acc}")
            else:
                with Entrez.esearch(db="nucleotide", term=f'"{config["organism"]}"[Organism] AND refseq[Filter]', retmax=5) as nh:
                    nuc_res = Entrez.read(nh)
                if nuc_res["IdList"]:
                    with Entrez.esummary(db="nucleotide", id=nuc_res["IdList"][0]) as nsh:
                        nuc_sum = Entrez.read(nsh)
                    found_acc = nuc_sum[0].get("Caption") or nuc_sum[0].get("AccessionVersion")
                    st.session_state["ref_acc_input"] = found_acc
                    st.sidebar.success(f"Found RefSeq Nucleotide: {found_acc}")
                else:
                    st.sidebar.warning("No RefSeq reference found automatically. Please enter an accession manually.")
        except Exception as e:
            st.sidebar.error(f"Lookup failed: {e}")

ref_acc = st.sidebar.text_input("Reference Accession", key="ref_acc_input")
config["reference_accession"] = ref_acc

if data_source == "NCBI":
    new_email = st.sidebar.text_input("NCBI Email", config.get("ncbi_email", ""))

st.sidebar.markdown("---")
new_max = st.sidebar.slider("Max Sequences to Download", 1, 50, config.get("max_sequences", 15))

# --- ADVANCED PARAMETERS ---
st.sidebar.subheader("🔬 Advanced Parameters")

new_min_len = st.sidebar.number_input(
    "Minimum Length (bp)", 
    value=int(config.get("min_length", auto_len))
)

new_max_n = st.sidebar.slider(
    "Max Ambiguous Bases (N%)", 
    0.0, 10.0, 
    float(config.get("max_ambiguous_n_percent", 1.0)), 
    0.1
)

align_options = {
    "Very Strict (<5% divergence)": "asm5",
    "Moderate (<10% divergence)": "asm10",
    "Relaxed (<20% divergence)": "asm20",
    "Highly Divergent / Short Segments (map-ont)": "map-ont"
}
current_preset = config.get("alignment_preset", "asm5")
current_label = [k for k, v in align_options.items() if v == current_preset]
default_idx = list(align_options.keys()).index(current_label[0]) if current_label else 0

align_choice = st.sidebar.selectbox("Alignment Stringency", list(align_options.keys()), index=default_idx)
new_align_preset = align_options[align_choice]

# 3 --- RUN & FLUSH BUTTONS ---
run_col, flush_col = st.sidebar.columns(2)

if run_col.button("🚀 Run Pipeline", width="stretch"):
    config["reference_accession"] = ref_acc
    config["search_mode"] = search_mode
    config["search_term"] = new_term
    config["ena_query"] = ena_query
    config["bvbrc_query"] = bvbrc_query   # bv-brc query
    config["ncbi_email"] = new_email
    config["max_sequences"] = new_max
    config["min_length"] = new_min_len
    config["max_ambiguous_n_percent"] = new_max_n
    config["alignment_preset"] = new_align_preset

    with open("config/config.yaml", "w") as f:
        yaml.dump(config, f)

    with st.spinner("Running Iostrack Engine..."):
        result = subprocess.run(
            [
                "snakemake",
                "-s",
                "workflow/Snakefile",
                "-c",
                "4",
                "--use-conda",
                "--forceall",
            ],
            capture_output=True,
            text=True,
        )

        if result.returncode == 0:
            st.sidebar.success("✅ Pipeline finished successfully!")
        else:
            st.sidebar.error("❌ Pipeline failed. Check terminal for details.")
            st.sidebar.text(result.stderr)

with flush_col.popover("🗑️ Flush", use_container_width=True):
    st.markdown("Are you sure you want to delete all files in `results/`?")
    if st.button("⚠️ Yes, Delete All", type="primary", width="stretch"):
        if os.path.exists("results"):
            shutil.rmtree("results")
            st.rerun()

# 4. The Main Results Dashboard
st.header("📊 Variant & Structural Analysis")

results_file = "results/summary/sample_report.csv"
splits_file = "results/summary/split_alignments.csv"
hotspot_file = "results/summary/snp_coordinates.csv"
genes_file = "results/summary/gene_features.csv"
vcf_file = "results/summary/variants.vcf"

if os.path.exists(results_file):
    df = pd.read_csv(results_file)
    
    if not df.empty:
        df = df.dropna(subset=["SNPs"])
        df = df.rename(columns={
            "Length_bp_x": "Length (bp)",
            "Length_bp": "Length (bp)"
        })
        
        # Load SNP coordinate data and GenBank gene features once
        snp_full_df = pd.DataFrame()
        if os.path.exists(hotspot_file):
            try:
                snp_full_df = pd.read_csv(hotspot_file)
                if not snp_full_df.empty:
                    snp_full_df["Accession"] = snp_full_df["Accession"].astype(str).str.strip()
            except pd.errors.EmptyDataError:
                snp_full_df = pd.DataFrame()
                
        genes_df = pd.DataFrame()
        if os.path.exists(genes_file):
            try:
                genes_df = pd.read_csv(genes_file)
            except pd.errors.EmptyDataError:
                genes_df = pd.DataFrame()
        
        # Calculate Cohort-Wide Ti/Tv and Missense/Synonymous Ratios
        total_ti = df["Transitions"].sum() if "Transitions" in df.columns else 0
        total_tv = df["Transversions"].sum() if "Transversions" in df.columns else 0
        cohort_titv = round(total_ti / total_tv, 2) if total_tv > 0 else float(total_ti)
        
        total_mis = df["Missense_SNPs"].sum() if "Missense_SNPs" in df.columns else 0
        total_syn = df["Synonymous_SNPs"].sum() if "Synonymous_SNPs" in df.columns else 0
        cohort_dnds = round(total_mis / total_syn, 2) if total_syn > 0 else float(total_mis)
            
        col1, col2, col3, col4, col5 = st.columns(5)
        col1.metric("Genomes Mapped", len(df))
        col2.metric("Average SNPs", round(df["SNPs"].mean(), 1) if not df.empty else 0)
        col3.metric("Max SNPs in a Sample", int(df["SNPs"].max()) if not df.empty else 0)
        col4.metric(
            "Cohort Ti/Tv Ratio",
            cohort_titv,
            help="Ratio of Transitions (A↔G, C↔T) to Transversions. Viral genomes typically show a strong transition bias (>2.0)."
        )
        col5.metric(
            "Missense / Syn Ratio",
            cohort_dnds,
            help=f"Ratio of amino-acid-altering SNPs ({total_mis}) to silent coding SNPs ({total_syn}) across the cohort."
        )
        
        tab1, tab2, tab3 = st.tabs([
            "🧬 Whole Genome Overview",
            "✂️ Structural Segments & Splits",
            "🌳 Clade Clustering & Phylogenomics"
        ])
        
        current_source = config.get("data_source", "NCBI")
        def make_url(acc):
            if current_source == "ENA":
                return f"https://www.ebi.ac.uk/ena/browser/view/{acc}"
            elif current_source == "BV-BRC":
                # BV-BRC uses GenBank accessions natively, so we route to NCBI
                # to avoid BV-BRC's heavy frontend failing to render.
                return f"https://www.ncbi.nlm.nih.gov/nuccore/{acc}"
            else:
                return f"https://www.ncbi.nlm.nih.gov/nuccore/{acc}"
        
        with tab1:
            all_accessions = df["Accession"].astype(str).unique().tolist()
            selected_isolates = st.multiselect(
                "🔍 Filter Substitution Spectrum & Highlights by Isolate Accession (leave empty for Entire Cohort):",
                options=all_accessions,
                default=[]
            )
            
            filtered_snps = (
                snp_full_df[snp_full_df["Accession"].isin(selected_isolates)]
                if (selected_isolates and not snp_full_df.empty) else snp_full_df
            )
            
            chart_col1, chart_col2 = st.columns([3, 2])
            
            with chart_col1:
                st.subheader("Total Mutational Burden")
                plot_df = df[df["Accession"].isin(selected_isolates)] if selected_isolates else df
                hover_cols = [
                    c for c in [
                        "UpdateDate", "Length (bp)", "Transitions", "Transversions", "Ti_Tv_Ratio",
                        "Missense_SNPs", "Synonymous_SNPs", "Missense_Syn_Ratio", "Insertions", "Deletions"
                    ] if c in plot_df.columns
                ]
                
                fig = px.bar(
                    plot_df, x="Accession", y="SNPs", color="GC_Percent",
                    hover_data=hover_cols,
                    title="Aggregated SNPs relative to Reference Genome"
                )
                st.plotly_chart(fig, width="stretch")
                
            with chart_col2:
                st.subheader("Substitution Spectrum (Ti/Tv)")
                if not filtered_snps.empty and "Mutation" in filtered_snps.columns:
                    sel_ti = (filtered_snps["Mut_Class"] == "Transition").sum()
                    sel_tv = (filtered_snps["Mut_Class"] == "Transversion").sum()
                    sel_ratio = round(sel_ti / sel_tv, 2) if sel_tv > 0 else float(sel_ti)
                    scope_label = f"{len(selected_isolates)} Selected Isolate(s)" if selected_isolates else "All Isolates"
                    
                    spec_counts = (
                        filtered_snps.groupby(["Mutation", "Mut_Class"])
                        .size()
                        .reset_index(name="Count")
                        .sort_values("Count", ascending=False)
                    )
                    fig_spec = px.bar(
                        spec_counts,
                        x="Mutation",
                        y="Count",
                        color="Mut_Class",
                        color_discrete_map={"Transition": "#00B4D8", "Transversion": "#F72585"},
                        title=f"{scope_label} | Ti/Tv = {sel_ratio} ({sel_ti} Ti : {sel_tv} Tv)"
                    )
                    fig_spec.update_layout(xaxis_title="Base Substitution (Ref > Alt)", yaxis_title="Total Occurrences")
                    st.plotly_chart(fig_spec, width="stretch")
                else:
                    st.info("Re-run pipeline to populate substitution spectrum.")
                    
            # --- Functional Consequence & Recurrent Amino Acid Mutations Panel ---
            if not filtered_snps.empty and "Consequence" in filtered_snps.columns:
                func_col1, func_col2 = st.columns([2, 3])
                with func_col1:
                    st.subheader("Functional Consequences")
                    conseq_counts = (
                        filtered_snps.groupby(["Gene", "Consequence"])
                        .size()
                        .reset_index(name="SNPs")
                    )
                    fig_conseq = px.bar(
                        conseq_counts,
                        x="Gene",
                        y="SNPs",
                        color="Consequence",
                        color_discrete_map={
                            "Synonymous": "#2EC4B6",
                            "Missense": "#FF9F1C",
                            "Nonsense (Stop)": "#E71D36",
                            "Non-coding": "#8D99AE"
                        },
                        title="Coding vs. Non-Coding Burden by Viral Gene"
                    )
                    fig_conseq.update_layout(xaxis_title="Reference Viral Gene", yaxis_title="SNP Count", height=340)
                    st.plotly_chart(fig_conseq, width="stretch")
                    
                with func_col2:
                    st.subheader("Top Recurrent Amino Acid Substitutions")
                    missense_only = filtered_snps[filtered_snps["Consequence"].isin(["Missense", "Nonsense (Stop)"])]
                    if not missense_only.empty:
                        top_aa = (
                            missense_only.groupby(["AA_Change", "Gene", "Product", "Mutation", "Consequence"])["Accession"]
                            .nunique()
                            .reset_index(name="Isolate_Count")
                            .sort_values("Isolate_Count", ascending=False)
                            .head(15)
                        )
                        top_aa["Cohort_Freq_%"] = ((top_aa["Isolate_Count"] / max(1, len(plot_df))) * 100).round(1)
                        st.dataframe(top_aa, width="stretch", height=300)
                    else:
                        st.info("No missense or nonsense amino acid changes detected in selection.")
            
            disp_df = plot_df.copy()
            disp_df["Accession"] = disp_df["Accession"].apply(make_url)
            
            st.dataframe(
                disp_df, 
                width="stretch",
                column_config={"Accession": st.column_config.LinkColumn("Accession", display_text=r"([A-Z]{1,6}_?[0-9]{4,})")}
            )
            
            btn_col1, btn_col2, _ = st.columns([1, 1, 3])
            with btn_col1:
                st.download_button(
                    "⬇️ Download Aggregated CSV",
                    df.to_csv(index=False).encode("utf-8"),
                    "iostrack_aggregated.csv",
                    "text/csv",
                    width="stretch"
                )
            with btn_col2:
                if os.path.exists(vcf_file):
                    with open(vcf_file, "rb") as vf:
                        st.download_button(
                            "⬇️ Download Multi-Sample VCF",
                            vf.read(),
                            "iostrack_variants.vcf",
                            "text/vcf",
                            width="stretch"
                        )
            
        with tab2:
            if os.path.exists(splits_file):
                splits_df = pd.read_csv(splits_file)
                if not splits_df.empty:
                    splits_df["Accession"] = splits_df["Accession"].astype(str).str.strip()
                
                # --- 1. Genomic Hotspot Map with Multi-Contig & Non-Overlapping CDS Support ---
                if not snp_full_df.empty:
                    snp_df = snp_full_df.copy()
                    st.subheader("🧬 Genomic Hotspot Map")
                    
                    has_contigs = "Target_Contig" in snp_df.columns and "Contig" in genes_df.columns
                    contig_list = sorted(snp_df["Target_Contig"].dropna().unique().tolist()) if has_contigs else []
                    
                    ctrl_col1, ctrl_col2, ctrl_col3 = st.columns([2, 2, 1])
                    with ctrl_col1:
                        view_mode = st.radio(
                            "Visualization Style",
                            ["🔥 2D Heatmap (Best for 10+ genomes)", "📊 Faceted Bar Chart (Detailed per-isolate view)"],
                            horizontal=True
                        )
                    with ctrl_col2:
                        if len(contig_list) > 1:
                            selected_contig = st.selectbox(
                                "Reference Contig / Segment",
                                ["All Contigs (Concatenated View)"] + contig_list
                            )
                        else:
                            selected_contig = contig_list[0] if contig_list else "All"
                    with ctrl_col3:
                        show_genes = st.checkbox("🧬 Overlay CDS", value=True)
                    
                    # Build coordinate space: either single contig or concatenated multi-contig offsets
                    plot_genes_df = genes_df.copy()
                    plot_splits_df = splits_df.copy()
                    
                    if has_contigs and len(contig_list) > 1:
                        if selected_contig != "All Contigs (Concatenated View)":
                            snp_df = snp_df[snp_df["Target_Contig"] == selected_contig].copy()
                            plot_genes_df = plot_genes_df[plot_genes_df["Contig"] == selected_contig].copy()
                            if "Target_Contig" in plot_splits_df.columns:
                                plot_splits_df = plot_splits_df[plot_splits_df["Target_Contig"] == selected_contig].copy()
                            snp_df["Plot_Pos"] = snp_df["Target_Pos"]
                        else:
                            # Calculate cumulative offsets so each segment sits side-by-side without overlapping
                            all_ctgs = sorted(
                                set(contig_list)
                                | (set(plot_genes_df["Contig"].dropna().unique()) if not plot_genes_df.empty else set())
                            )
                            ctg_offsets = {}
                            current_offset = 0
                            gap = 150  # Visual spacer between segments
                            for ctg in all_ctgs:
                                ctg_offsets[ctg] = current_offset
                                max_snp = snp_df[snp_df["Target_Contig"] == ctg]["Target_Pos"].max() if not snp_df[snp_df["Target_Contig"] == ctg].empty else 0
                                max_gene = plot_genes_df[plot_genes_df["Contig"] == ctg]["End"].max() if not plot_genes_df[plot_genes_df["Contig"] == ctg].empty else 0
                                max_split = (
                                    plot_splits_df[plot_splits_df["Target_Contig"] == ctg]["Target_End"].max()
                                    if ("Target_Contig" in plot_splits_df.columns and not plot_splits_df[plot_splits_df["Target_Contig"] == ctg].empty)
                                    else 0
                                )
                                ctg_len = int(max(max_snp + 1, max_gene, max_split, 100))
                                current_offset += ctg_len + gap
                                
                            snp_df["Plot_Pos"] = snp_df.apply(
                                lambda r: int(r["Target_Pos"]) + ctg_offsets.get(r["Target_Contig"], 0), axis=1
                            )
                            if not plot_genes_df.empty:
                                plot_genes_df["Start"] = plot_genes_df.apply(
                                    lambda r: int(r["Start"]) + ctg_offsets.get(r["Contig"], 0), axis=1
                                )
                                plot_genes_df["End"] = plot_genes_df.apply(
                                    lambda r: int(r["End"]) + ctg_offsets.get(r["Contig"], 0), axis=1
                                )
                    else:
                        snp_df["Plot_Pos"] = snp_df["Target_Pos"]
                    
                    max_genome_pos = int(max(
                        snp_df["Plot_Pos"].max() + 1 if not snp_df.empty else 100,
                        plot_genes_df["End"].max() if not plot_genes_df.empty else 0,
                        plot_splits_df["Target_End"].max() if (not plot_splits_df.empty and selected_contig != "All Contigs (Concatenated View)") else 0
                    ))
                    num_bins = st.slider(
                        "🔍 Heatmap Resolution (Number of Genomic Windows)",
                        min_value=30, max_value=300, value=120, step=30,
                        help="Increase to 180–300 windows to separate small adjacent genes like the 228-bp SARS-CoV-2 E gene."
                    )
                    bin_edges = [int(i * max_genome_pos / num_bins) for i in range(num_bins + 1)]
                    bin_edges = [max(bin_edges[i], i) for i in range(num_bins + 1)]
                    all_window_mids = [int((bin_edges[i] + bin_edges[i + 1]) / 2) for i in range(num_bins)]
                    
                    bins = pd.cut(snp_df["Plot_Pos"], bins=bin_edges, labels=all_window_mids, include_lowest=True)
                    snp_df["Window"] = bins.astype(float).fillna(0).astype(int)
                    
                    binned_df = snp_df.groupby(["Accession", "Window"]).size().reset_index(name="SNPs")
                    num_accessions = max(1, len(binned_df["Accession"].unique()))
                    
                    x_label = (
                        "Concatenated Genome Position Across Segments (bp)"
                        if selected_contig == "All Contigs (Concatenated View)"
                        else f"Reference Position on {selected_contig} (bp)"
                    )
                    
                    if "Heatmap" in view_mode:
                        st.markdown("Each row is an isolate and each column is a genomic window. **Darker red cells** reveal shared mutation hotspots across the cohort.")
                        
                        heatmap_matrix = (
                            binned_df.pivot(index="Accession", columns="Window", values="SNPs")
                            .reindex(columns=all_window_mids, fill_value=0)
                            .fillna(0)
                        )
                        
                        fig_map = px.imshow(
                            heatmap_matrix,
                            labels=dict(x=x_label, y="Isolate Accession", color="SNPs"),
                            color_continuous_scale="Reds",
                            aspect="auto",
                            title="SNP Density & Frequency Across Reference Genome"
                        )
                        fig_map.update_layout(height=max(450, 180 + (num_accessions * 24)))
                        fig_map.update_xaxes(tickformat="d")
                    else:
                        st.markdown("This chart divides the reference genome into windows. **Bar height** and **darker red colors** indicate hyper-mutated genomic regions for each isolate.")
                        
                        safe_spacing = min(0.03, 0.5 / max(1, num_accessions - 1))
                        
                        fig_map = px.bar(
                            binned_df,
                            x="Window",
                            y="SNPs",
                            color="SNPs",
                            facet_row="Accession",
                            facet_row_spacing=safe_spacing,
                            color_continuous_scale="Reds",
                            title="SNP Density & Frequency Across Reference Genome"
                        )
                        fig_map.update_layout(
                            xaxis_title=x_label,
                            height=280 + (num_accessions * 80)
                        )
                        fig_map.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
                        fig_map.update_xaxes(tickformat="d", range=[0, max_genome_pos])
                        
                    # Overlay GenBank CDS boundaries with automatic overlap merging
                    if show_genes and not plot_genes_df.empty:
                        group_cols = ["Contig", "Gene"] if "Contig" in plot_genes_df.columns else ["Gene"]
                        raw_genes = (
                            plot_genes_df.groupby(group_cols, as_index=False, sort=False)
                            .agg({"Start": "min", "End": "max", "Product": "first"})
                            .sort_values(["Start", "End"], ascending=[True, False])
                        )
                        
                        # Merge nested/overlapping ORFs (e.g. M1/M2, PB1/PB1-F2) so vrects never overlap
                        merged_intervals = []
                        for _, row in raw_genes.iterrows():
                            s, e, g = int(row["Start"]), int(row["End"]), str(row["Gene"])
                            if not merged_intervals:
                                merged_intervals.append({"Start": s, "End": e, "Genes": [g]})
                            else:
                                prev = merged_intervals[-1]
                                overlap = min(prev["End"], e) - max(prev["Start"], s)
                                shorter_len = max(1, min(prev["End"] - prev["Start"], e - s))
                                
                                if overlap > 0:
                                    # If heavily overlapping/nested (>=20%), combine into one clean block
                                    if (overlap / shorter_len) >= 0.20:
                                        prev["End"] = max(prev["End"], e)
                                        if g not in prev["Genes"]:
                                            prev["Genes"].append(g)
                                    else:
                                        # Minor terminal overlap: split cleanly at the midpoint
                                        mid_split = int((max(prev["Start"], s) + min(prev["End"], e)) / 2)
                                        prev["End"] = mid_split
                                        merged_intervals.append({"Start": mid_split + 1, "End": e, "Genes": [g]})
                                else:
                                    merged_intervals.append({"Start": s, "End": e, "Genes": [g]})
                        
                        fill_palette = ["rgba(0, 180, 216, 0.14)", "rgba(247, 37, 133, 0.10)"]
                        line_palette = ["rgba(0, 180, 216, 0.85)", "rgba(247, 37, 133, 0.85)"]
                        badge_palette = ["#0077B6", "#C9184A"]
                        y_shifts = [2, 18, 34]
                        
                        for idx, g_row in enumerate(merged_intervals):
                            g_start = g_row["Start"]
                            g_end = g_row["End"]
                            g_mid = int((g_start + g_end) / 2)
                            label_text = "/".join(g_row["Genes"])
                            
                            fig_map.add_vrect(
                                x0=g_start,
                                x1=g_end,
                                fillcolor=fill_palette[idx % 2],
                                layer="above",
                                line_width=1.5,
                                line_dash="dot",
                                line_color=line_palette[idx % 2]
                            )
                            
                            fig_map.add_annotation(
                                x=g_mid,
                                y=1.0,
                                xref="x",
                                yref="paper",
                                text=f"<b>{label_text}</b>",
                                showarrow=False,
                                yshift=y_shifts[idx % 3],
                                font=dict(size=11, color="#FFFFFF"),
                                bgcolor=badge_palette[idx % 2],
                                borderpad=3,
                                opacity=0.92
                            )
                            
                        fig_map.update_layout(margin=dict(t=95))
                    
                    st.plotly_chart(fig_map, width="stretch")
                    
                    if show_genes and not genes_df.empty:
                        with st.expander("📖 View Reference Viral Gene (CDS) Coordinates Table"):
                            st.dataframe(genes_df, width="stretch")
                            
                    st.markdown("---")
                    
                    # --- 2. Segment, Gene & Genomic Window Ti/Tv Inspector ---
                    if "Mutation" in snp_full_df.columns and not splits_df.empty:
                        st.subheader("🔬 Segment, Viral Gene & Window Ti/Tv Inspector")
                        st.markdown("Inspect the transition/transversion signature and coding impact inside a specific **viral gene (CDS)**, **aligned segment (chunk)**, or custom **coordinate window**.")
                        
                        insp_col1, insp_col2 = st.columns([1, 2])
                        
                        with insp_col1:
                            filter_modes = ["By Alignment Segment (Chunk)", "By Genomic Coordinate Window (bp)"]
                            if not genes_df.empty:
                                filter_modes.insert(0, "By Viral Gene (CDS)")
                                
                            chunk_mode = st.radio("Filter Mode", filter_modes)
                            
                            if chunk_mode == "By Viral Gene (CDS)":
                                if "Contig" in genes_df.columns:
                                    genes_df["Gene_Label"] = (
                                        genes_df["Contig"].astype(str) + " | " + genes_df["Gene"].astype(str) +
                                        " (" + genes_df["Product"].astype(str) + " | " +
                                        genes_df["Start"].astype(str) + "–" + genes_df["End"].astype(str) + " bp)"
                                    )
                                else:
                                    genes_df["Gene_Label"] = (
                                        genes_df["Gene"].astype(str) + " (" + genes_df["Product"].astype(str) +
                                        " | " + genes_df["Start"].astype(str) + "–" + genes_df["End"].astype(str) + " bp)"
                                    )
                                chosen_gene_label = st.selectbox("Select Reference Viral Gene:", genes_df["Gene_Label"].tolist())
                                gene_row = genes_df[genes_df["Gene_Label"] == chosen_gene_label].iloc[0]
                                
                                gene_isolates = st.multiselect(
                                    "Filter Gene by Isolate (optional):",
                                    options=splits_df["Accession"].unique().tolist(),
                                    default=[]
                                )
                                chunk_snps = snp_full_df[
                                    (snp_full_df["Target_Pos"] + 1 >= int(gene_row["Start"])) &
                                    (snp_full_df["Target_Pos"] + 1 <= int(gene_row["End"]))
                                ]
                                if "Contig" in gene_row and "Target_Contig" in chunk_snps.columns:
                                    chunk_snps = chunk_snps[chunk_snps["Target_Contig"] == gene_row["Contig"]]
                                if gene_isolates:
                                    chunk_snps = chunk_snps[chunk_snps["Accession"].isin(gene_isolates)]
                                chunk_title = f"Gene {gene_row['Gene']} ({gene_row['Start']}–{gene_row['End']} bp)"
                                
                            elif chunk_mode == "By Alignment Segment (Chunk)":
                                contig_prefix = (splits_df["Target_Contig"].astype(str) + " | ") if "Target_Contig" in splits_df.columns else ""
                                splits_df["Segment_Label"] = (
                                    splits_df["Accession"] + " | " + contig_prefix + "Ref: " +
                                    splits_df["Target_Start"].astype(str) + "–" +
                                    splits_df["Target_End"].astype(str) + " bp (" +
                                    splits_df["Strand"].astype(str) + ")"
                                )
                                seg_options = ["All Aligned Segments"] + splits_df["Segment_Label"].tolist()
                                chosen_seg = st.selectbox("Select Alignment Chunk:", seg_options)
                                
                                if chosen_seg == "All Aligned Segments":
                                    chunk_snps = snp_full_df
                                    chunk_title = "All Segments"
                                else:
                                    seg_row = splits_df[splits_df["Segment_Label"] == chosen_seg].iloc[0]
                                    chunk_snps = snp_full_df[
                                        (snp_full_df["Accession"] == seg_row["Accession"]) &
                                        (snp_full_df["Target_Pos"] >= seg_row["Target_Start"]) &
                                        (snp_full_df["Target_Pos"] <= seg_row["Target_End"])
                                    ]
                                    if "Target_Contig" in seg_row and "Target_Contig" in chunk_snps.columns:
                                        chunk_snps = chunk_snps[chunk_snps["Target_Contig"] == seg_row["Target_Contig"]]
                                    chunk_title = chosen_seg
                            else:
                                min_pos = int(splits_df["Target_Start"].min())
                                max_pos = int(splits_df["Target_End"].max())
                                coord_range = st.slider(
                                    "Reference Coordinate Window (bp):",
                                    min_value=min_pos,
                                    max_value=max(min_pos + 1, max_pos),
                                    value=(min_pos, max_pos),
                                    step=100
                                )
                                chunk_isolates = st.multiselect(
                                    "Filter Window by Isolate (optional):",
                                    options=splits_df["Accession"].unique().tolist(),
                                    default=[]
                                )
                                chunk_snps = snp_full_df[
                                    (snp_full_df["Target_Pos"] >= coord_range[0]) &
                                    (snp_full_df["Target_Pos"] <= coord_range[1])
                                ]
                                if chunk_isolates:
                                    chunk_snps = chunk_snps[chunk_snps["Accession"].isin(chunk_isolates)]
                                chunk_title = f"Window {coord_range[0]:,}–{coord_range[1]:,} bp"
                                
                            c_ti = int((chunk_snps["Mut_Class"] == "Transition").sum()) if not chunk_snps.empty else 0
                            c_tv = int((chunk_snps["Mut_Class"] == "Transversion").sum()) if not chunk_snps.empty else 0
                            c_ratio = round(c_ti / c_tv, 2) if c_tv > 0 else float(c_ti)
                            
                            m_col1, m_col2, m_col3 = st.columns(3)
                            m_col1.metric("Selection Ti/Tv", c_ratio)
                            m_col2.metric("Transitions", c_ti)
                            m_col3.metric("Transversions", c_tv)
                            
                            if "Consequence" in chunk_snps.columns:
                                c_mis = int(chunk_snps["Consequence"].isin(["Missense", "Nonsense (Stop)"]).sum())
                                c_syn = int((chunk_snps["Consequence"] == "Synonymous").sum())
                                c_dnds = round(c_mis / c_syn, 2) if c_syn > 0 else float(c_mis)
                                k_col1, k_col2, k_col3 = st.columns(3)
                                k_col1.metric("Missense/Syn", c_dnds)
                                k_col2.metric("Missense", c_mis)
                                k_col3.metric("Synonymous", c_syn)
                            
                        with insp_col2:
                            if not chunk_snps.empty:
                                chunk_spec = (
                                    chunk_snps.groupby(["Mutation", "Mut_Class"])
                                    .size()
                                    .reset_index(name="Count")
                                    .sort_values("Count", ascending=False)
                                )
                                fig_chunk_spec = px.bar(
                                    chunk_spec,
                                    x="Mutation",
                                    y="Count",
                                    color="Mut_Class",
                                    color_discrete_map={"Transition": "#00B4D8", "Transversion": "#F72585"},
                                    title=f"Substitution Spectrum: {chunk_title}"
                                )
                                fig_chunk_spec.update_layout(
                                    xaxis_title="Base Substitution (Ref > Alt)",
                                    yaxis_title="Occurrences in Selection",
                                    height=340
                                )
                                st.plotly_chart(fig_chunk_spec, width="stretch")
                            else:
                                st.info("No SNPs found inside the selected gene/chunk/window.")
                        with st.expander("🧬 Single-Nucleotide & Codon-Level Mutation Table (1-Based Coordinates)", expanded=False):
                            if not chunk_snps.empty:
                                single_nt_df = chunk_snps.copy()
                                single_nt_df["Position_1Based"] = single_nt_df["Target_Pos"] + 1
                                nt_cols = [
                                    c for c in [
                                        "Accession", "Target_Contig", "Position_1Based",
                                        "Mutation", "Mut_Class", "Gene", "Consequence",
                                        "Codon_Change", "AA_Change", "Product"
                                    ] if c in single_nt_df.columns
                                ]
                                st.dataframe(
                                    single_nt_df[nt_cols].sort_values(["Position_1Based", "Accession"]),
                                    width="stretch",
                                    height=320
                                )
                            else:
                                st.info("No SNPs in current selection.")
                        st.markdown("---")
                
                # --- 3. Alignment Segments Table ---
                st.subheader("✂️ Alignment Segments & Coordinate Data")
                st.markdown("This table shows the exact Reference `Target` and Isolate `Query` start/end coordinates for every mapped segment. Secondary overlapping alignments have been filtered out.")
                
                disp_splits = splits_df.copy()
                disp_splits["Accession"] = disp_splits["Accession"].apply(make_url)
                
                cols_to_show = [
                    c for c in [
                        "Accession", "Target_Contig", "Target_Start", "Target_End",
                        "Query_Start", "Query_End", "Strand",
                        "SNPs", "Transitions", "Transversions", "Ti_Tv_Ratio",
                        "Missense_SNPs", "Synonymous_SNPs", "Noncoding_SNPs",
                        "Insertions", "Deletions"
                    ] if c in disp_splits.columns
                ]
                st.dataframe(
                    disp_splits[cols_to_show], 
                    width="stretch",
                    column_config={"Accession": st.column_config.LinkColumn("Accession", display_text=r"([A-Z]{1,6}_?[0-9]{4,})")}
                )
                st.download_button(
                    "⬇️ Download Coordinate CSV",
                    splits_df.drop(columns=["Segment_Label"], errors="ignore").to_csv(index=False).encode('utf-8'),
                    "iostrack_coordinates.csv",
                    "text/csv"
                )

        with tab3:
            st.subheader("🌳 Pairwise SNP Distance Matrix & Clade Clustering")
            st.markdown(
                "This matrix computes the exact **pairwise nucleotide substitution distance** between every mapped isolate "
                "and orders samples by genetic similarity so viral lineages/subgroups cluster into distinct blocks."
            )
            
            if not snp_full_df.empty and len(df) >= 2:
                clust_col1, clust_col2 = st.columns([1, 2])
                with clust_col1:
                    include_ref = st.checkbox("Include Reference Genome in Matrix", value=True)
                with clust_col2:
                    max_dist_thresh = st.slider(
                        "🔍 Outbreak / Lineage Link Threshold (Max Pairwise SNPs):",
                        min_value=0,
                        max_value=int(max(50, df["SNPs"].max() * 2)),
                        value=min(25, int(max(10, df["SNPs"].median()))),
                        help="Isolate pairs differing by fewer than this many SNPs are listed in the Linked Pairs table below."
                    )
                
                # Build exact allele sets per isolate: {((Target_Contig, Target_Pos), Alt_Base), ...}
                sample_names = df["Accession"].astype(str).tolist()
                if include_ref:
                    ref_label = f"REF ({config.get('reference_accession', 'Ref')})"
                    sample_names = [ref_label] + sample_names
                else:
                    ref_label = None
                    
                snp_sets = {s: set() for s in sample_names}
                for _, s_row in snp_full_df.iterrows():
                    acc = str(s_row["Accession"])
                    ctg = str(s_row.get("Target_Contig", "Ref"))
                    if acc in snp_sets:
                        snp_sets[acc].add(((ctg, int(s_row["Target_Pos"])), str(s_row["Alt_Base"])))
                        
                n_samp = len(sample_names)
                dist_matrix = pd.DataFrame(0, index=sample_names, columns=sample_names)
                
                for i in range(n_samp):
                    for j in range(i + 1, n_samp):
                        s1, s2 = sample_names[i], sample_names[j]
                        diff_loci = snp_sets[s1].symmetric_difference(snp_sets[s2])
                        pos_count = len({locus for locus, _ in diff_loci})
                        dist_matrix.loc[s1, s2] = pos_count
                        dist_matrix.loc[s2, s1] = pos_count
                        
                ordered_samples = [sample_names[0]]
                remaining = set(sample_names[1:])
                while remaining:
                    last = ordered_samples[-1]
                    next_s = min(remaining, key=lambda x: dist_matrix.loc[last, x])
                    ordered_samples.append(next_s)
                    remaining.remove(next_s)
                    
                dist_matrix = dist_matrix.loc[ordered_samples, ordered_samples]
                
                fig_dist = px.imshow(
                    dist_matrix,
                    labels=dict(x="Isolate", y="Isolate", color="Pairwise SNPs"),
                    color_continuous_scale="Viridis",
                    aspect="auto",
                    title="Ordered Pairwise SNP Distance Matrix (Darker = More Closely Related)"
                )
                fig_dist.update_layout(height=max(480, 180 + (n_samp * 22)))
                st.plotly_chart(fig_dist, width="stretch")
                
                st.markdown("---")
                pair_col1, pair_col2 = st.columns([1, 1])
                
                with pair_col1:
                    st.subheader(f"🔗 Closely Linked Isolate Pairs (≤ {max_dist_thresh} SNPs)")
                    pairs_list = []
                    for i in range(len(ordered_samples)):
                        for j in range(i + 1, len(ordered_samples)):
                            s1, s2 = ordered_samples[i], ordered_samples[j]
                            if ref_label and (s1 == ref_label or s2 == ref_label):
                                continue
                            d = int(dist_matrix.loc[s1, s2])
                            if d <= max_dist_thresh:
                                shared = len(snp_sets[s1].intersection(snp_sets[s2]))
                                pairs_list.append({
                                    "Isolate_1": s1,
                                    "Isolate_2": s2,
                                    "Pairwise_SNP_Distance": d,
                                    "Shared_SNPs": shared
                                })
                    pairs_df = pd.DataFrame(pairs_list)
                    if not pairs_df.empty:
                        pairs_df = pairs_df.sort_values("Pairwise_SNP_Distance", ascending=True)
                        st.dataframe(pairs_df, width="stretch", height=300)
                    else:
                        st.info(f"No isolate pairs found within ≤ {max_dist_thresh} SNPs. Increase the threshold slider above.")
                        
                with pair_col2:
                    st.subheader("📌 Clade-Defining / High-Frequency SNPs")
                    loci_cols = [c for c in ["Target_Contig", "Target_Pos", "Mutation", "Gene", "Consequence", "AA_Change"] if c in snp_full_df.columns]
                    core_snps = (
                        snp_full_df.groupby(loci_cols)["Accession"]
                        .nunique()
                        .reset_index(name="Isolates_With_SNP")
                        .sort_values("Isolates_With_SNP", ascending=False)
                        .head(20)
                    )
                    core_snps["Target_Pos"] = core_snps["Target_Pos"] + 1
                    core_snps["Cohort_Freq_%"] = ((core_snps["Isolates_With_SNP"] / len(df)) * 100).round(1)
                    st.dataframe(core_snps, width="stretch", height=300)
                    
                st.download_button(
                    "⬇️ Download Pairwise Distance Matrix CSV",
                    dist_matrix.to_csv().encode("utf-8"),
                    "iostrack_pairwise_snp_distance.csv",
                    "text/csv"
                )
            else:
                st.info("At least 2 mapped genomes with SNP coordinates are required to compute pairwise clustering.")
    else:
        st.warning("The pipeline ran, but no sequences aligned successfully.")
else:
    st.info("No results found. Adjust parameters and click 'Run Pipeline'.")