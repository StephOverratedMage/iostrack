# 🧬 Iostrack: Automated Viral Genomics & Mutational Surveillance Engine

![Python](https://img.shields.io/badge/Python-3.11-blue?logo=python)
![Snakemake](https://img.shields.io/badge/Workflow-Snakemake-brightgreen?logo=snakemake)
![Streamlit](https://img.shields.io/badge/Dashboard-Streamlit-FF4B4B?logo=streamlit)
![Minimap2](https://img.shields.io/badge/Aligner-Minimap2-orange)
![Tests](https://img.shields.io/badge/Tests-Pytest-0A9EDC?logo=pytest)

**Iostrack** is an end-to-end bioinformatics workflow and interactive surveillance GUI designed to fetch, quality-filter, align, and functionally annotate viral isolate genomes from **NCBI Entrez** and the **European Nucleotide Archive (ENA)**.

Built with **Snakemake** for reproducible workflow orchestration and **Streamlit + Plotly** for real-time visual analytics, Iostrack transforms raw unaligned FASTA sequences into actionable epidemiological metrics—including **GenBank CDS amino acid consequences**, **$Ti/Tv$ mutational signatures**, **structural split-alignment coordinates**, **pairwise SNP distance matrices**, and **standards-compliant multi-sample VCFv4.2 exports**.

---

## 🏗️ Pipeline Architecture

```mermaid
flowchart LR
    UI["🖥️ Streamlit GUI\n(config.yaml)"] --> SM["🐍 Snakemake Router"]
    SM -->|data_source: NCBI| NCBI["fetch_ncbi.py\n(Entrez E-Utilities)"]
    SM -->|data_source: ENA| ENA["fetch_ena.py\n(ENA Portal & Browser API)"]
    
    NCBI --> RAW["results/raw/\nreference.gb & .fasta\nisolates.fasta\nmetadata.csv"]
    ENA --> RAW
    
    RAW --> QC["qc_filter.py\nLength & Ambiguous Base (N%) QC"]
    QC --> MM2["minimap2 --cs --secondary=no\nPairwise Reference Alignment"]
    
    MM2 --> PARSE["parse_mutations.py\ncs-Tag State Machine & CDS Translator"]
    RAW --> PARSE
    
    PARSE --> OUT1["📊 sample_report.csv\n(Burden, Ti/Tv, Missense/Syn)"]
    PARSE --> OUT2["🧬 snp_coordinates.csv\n(HGVS p. AA Changes & Loci)"]
    PARSE --> OUT3["✂️ split_alignments.csv\n(Structural Segments)"]
    PARSE --> OUT4["📄 variants.vcf\n(Multi-Sample VCFv4.2)"]