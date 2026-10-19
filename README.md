# 🧬 Iostrack: Automated Viral Genomics & Mutational Surveillance Engine

![Python](https://img.shields.io/badge/Python-3.11-blue?logo=python)
![Snakemake](https://img.shields.io/badge/Workflow-Snakemake-brightgreen?logo=snakemake)
![Streamlit](https://img.shields.io/badge/Dashboard-Streamlit-FF4B4B?logo=streamlit)
![Minimap2](https://img.shields.io/badge/Aligner-Minimap2-orange)
![Tests](https://img.shields.io/badge/Tests-Pytest-0A9EDC?logo=pytest)

**Iostrack** is an end-to-end bioinformatics workflow and interactive surveillance GUI designed to fetch, quality-filter, align, and functionally annotate viral isolate genomes from **NCBI Entrez** and the **European Nucleotide Archive (ENA)**.

Built with **Snakemake** for reproducible workflow orchestration and **Streamlit + Plotly** for real-time visual analytics, Iostrack transforms raw unaligned FASTA sequences into actionable epidemiological metrics—including **GenBank CDS amino acid consequences**, **$Ti/Tv$ mutational signatures**, **structural split-alignment coordinates**, **pairwise SNP distance matrices**, and **standards-compliant multi-sample VCFv4.2 exports**.

Iostrack is open source and free to use. Includes a custom Dark and Light mode depending on visual preference. Also supports WSL!

---
## 🖼️ GUI Environment

Interactive plots which visualize SNP counts and types across multiple sequences compared to the reference sequence.

<img width="930" height="793" alt="tab1" src="https://github.com/user-attachments/assets/5a99b093-3cb0-4dae-a09d-c9c7a465876b" />
<img width="942" height="805" alt="tab2" src="https://github.com/user-attachments/assets/d3d82461-c8c9-44b6-80aa-ee723271f9b7" />
<img width="966" height="822" alt="tab3" src="https://github.com/user-attachments/assets/d3947a3a-6be2-4479-b2db-eb9be2759066" />


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
```
---
## 📥 Installation

From your preferred Linux/WSL command prompt: 

```bash
git clone https://github.com/StephOverratedMage/iostrack.git
cd iostrack

conda env create -f environment.yml
conda activate iostrack

#Run the automated tests to verify installation
python -m pytest -v

#Open the interactive GUI dashboard
streamlit run app.py
```
#### Once the GUI opens, configure your settings in the sidebar to the left and then click "**🚀 Run Pipeline**" to start pipeline execution. NOTE: When querying NCBI, a valid e-mail address is required for their API.
---
## 🔜 Future prospects

**Iostrack is a prototype and work in progress. Features coming soon include: More customizability within plots, as well as more diverse and useful visualizations. More ways to analyze and sort results. 
Any feedback is appreciated!** 
