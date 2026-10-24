# 🧬 Iostrack: Automated Viral Genomics & Mutational Surveillance Engine

![Python](https://img.shields.io/badge/Python-3.11-blue?logo=python)
![Snakemake](https://img.shields.io/badge/Workflow-Snakemake-brightgreen?logo=snakemake)
![Streamlit](https://img.shields.io/badge/Dashboard-Streamlit-FF4B4B?logo=streamlit)
![Minimap2](https://img.shields.io/badge/Aligner-Minimap2-orange)
![Tests](https://img.shields.io/badge/Tests-Pytest-0A9EDC?logo=pytest)

**Iostrack** is an end-to-end bioinformatics workflow and interactive surveillance GUI for fetching, quality-filtering, aligning, and annotating viral isolate genomes. "Ios" comes from the Greek word meaning "virus".

It supports **NCBI Entrez** and the **European Nucleotide Archive (ENA)** as data sources.

Built with **Snakemake** for reproducible workflow orchestration and **Streamlit + Plotly** for real-time visual analytics, Iostrack transforms raw unaligned FASTA sequences into actionable epidemiological metrics—including **GenBank CDS amino acid consequences**, **$Ti/Tv$ mutational signatures**, **structural split-alignment coordinates**, **pairwise SNP distance matrices**, and **standards-compliant multi-sample VCFv4.2 exports**.

Iostrack is open source and free to use. It includes custom Dark and Light themes and supports Linux and WSL2 environments.

> **Status:** Research prototype. Suitable for exploratory analysis; not intended for clinical or diagnostic use.

---
## 📋 Contents
- [Definitions](#-definitions)
- [Features](#️-features)
- [GUI Environment](#️-gui-environment)
- [Pipeline Architecture](#️-pipeline-architecture)
- [Requirements](#-requirements)
- [Installation](#-installation)
- [Example Run](#-example-run)
- [Results](#-results)
- [Roadmap](#-roadmap)
- [License](#-license)

---
## 🔎 Definitions

* N content / N%: Ambiguous base call percentage.
* Ti/Tv: The ratio of transitions to transversions. Transitions are A↔G and C↔T substitutions; transversions are substitutions between a purine and a pyrimidine.
* CDS: Coding DNA Sequence. A segment of DNA that is transcribed and eventually translated into a protein molecule.
* Missense / Synonymous Ratio (dN/dS proxy): The ratio of amino-acid-altering (non-synonymous) SNPs to silent (synonymous) SNPs across coding regions. This is a simple mutation-count proxy and is not equivalent to a formal dN/dS estimate.
* PAF: Pairwise Alignment Format, a tab-delimited format for storing sequence alignment information.
* HGVS: Human Genome Variation Society nomenclature used for describing genetic variants.
* VCF: Variant Call Format, a standard text file format for storing DNA sequence variations.
* AF: Info field for VCF files, displaying allele frequency for each alternative allele.
* GT: Genotype

---
## 🛠️ Features

* Allows for querying the NCBI or ENA databases to fetch DNA sequences based on search terms (Organism, title terms) or accession numbers.
* Aligns all results to a reference genome specified by the user.
* Filters sequences by minimum length and maximum ambiguous-base percentage, provides configurable Minimap2 alignment presets for different expected divergence levels and limits the maximum number of sequences downloaded, displayed, and analyzed.
* Displays mutation summary tables and matrices detailing total mutations per accession or per primary segment based on the Minimap2 alignments. Mutation types (transition/transversion) as well as their ratio are also included.
* Visualizes total SNPs and GC content per accession.
* Produces a heatmap of mutations across the reference genome, with the option to overlay coding region limits over the heatmap.
* Visualizes the pairwise SNP distance between all analyzed genomes.
* Generates VCF, CSV and alignment reports.

## 🖼️ GUI Environment

Configure search and filtering through the sidebar:

<img width="330" height="598" alt="iostrack sidebar showing NCBI reference and isolate configuration" src="https://github.com/user-attachments/assets/948d42fa-847e-4f08-bd35-ccede671943b" />

The dashboard includes interactive plots that visualize SNP counts and mutation types across multiple sequences compared with the reference genome.

<img width="930" height="793" alt="iostrack mutation summary tab" src="https://github.com/user-attachments/assets/5a99b093-3cb0-4dae-a09d-c9c7a465876b" />
<img width="942" height="805" alt="iostrack mutation heatmap and smaller segment mutation table" src="https://github.com/user-attachments/assets/d3d82461-c8c9-44b6-80aa-ee723271f9b7" />
<img width="966" height="822" alt="iostrack pairwise SNP distance visualization" src="https://github.com/user-attachments/assets/d3947a3a-6be2-4479-b2db-eb9be2759066" />


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
## ✅ Requirements

* **Linux or WSL2** due to Minimap2 and Snakemake being intended for Linux. Native Windows is untested.
* Conda, Miniconda, Mamba or Miniforge. IDEs such as VSCode with WSL access also work.
* Python 3.11 (other versions not tested yet).
* Internet access for NCBI or ENA queries. If results are already present, it is not required.

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
## 💡 Example Run

1. Select `NCBI` as the data source.
2. Enter a reference accession, such as `NC_001781.1`.
3. Select **By Organism & Keyword**.
4. Enter the organism and sequence title filter. Leave the title filter empty to automatically disable it.
5. Provide an NCBI email address.
6. Set quality-control thresholds.
7. Click **🚀 Run Pipeline**.

> **Note:** NCBI Entrez queries require a valid email address for API identification and responsible usage.

---
## 📊 Results

Iostrack organizes its outputs hierarchically inside the `results/` directory. Below is the reference guide for every file generated during a pipeline execution:

| Directory / File | Category | Description & Downstream Purpose |
| :--- | :--- | :--- |
| **`results/raw/reference.fasta`** | Reference | Raw nucleotide sequence of the baseline reference genome used by `minimap2`. |
| **`results/raw/reference.gb`** | Reference | GenBank flat file (`.gb`) parsed by Biopython to extract viral gene boundaries (`CDS`) and reading frames. |
| **`results/raw/isolates.fasta`** | Raw Input | Pre-screened isolate genomes retained by the NCBI/ENA fetcher after length, ambiguous-base, duplicate, and reference-exclusion checks. |
| **`results/raw/metadata.csv`** | Metadata | Cohort metadata (`Accession`, `Title`, `Length_bp`, `UpdateDate`) mapped to GUI dataframes. |
| **`results/qc/clean_isolates.fasta`** | Quality Control | Final QC output passed to Minimap2 after the workflow’s QC rule. |
| **`results/qc/qc_report.csv`** | Quality Control | Per-sample QC metrics table tracking `Length_bp`, `GC_Percent`, `N_Percent`, and pass/fail `Status`. |
| **`results/align/alignments.paf`** | Alignment | Raw Pairwise Alignment Format (`.paf`) containing high-resolution mapping coordinates and `cs:Z:` tags. |
| **`results/summary/split_alignments.csv`** | Structural | Detailed coordinate breakdown of primary mapped segments. Secondary and supplementary alignment records are excluded. |
| **`results/summary/gene_features.csv`** | Annotation | Structured table of reference viral genes (`Gene`, `Product`, `Start`, `End`, `Strand`) for CDS track overlays. |
| **`results/summary/sample_report.csv`** | Summary Report | Aggregated cohort report combining mutational burdens, Ti/Tv ratios, and missense/synonymous mutation-count proxy metrics. |
| **`results/summary/snp_coordinates.csv`** | Variant Loci | Granular locus-by-locus table of every SNP, tracking base changes, mutational classes, and HGVS protein changes. |
| **`results/summary/variants.vcf`** | Interoperable | **Multi-Sample VCF v4.2** file with `NS`, `AF`, `GENE`, `CONSEQ`, `AA`, and binary sample genotypes (`GT`). |

> **Note on Re-running Queries:** Clicking **🚀 Run Pipeline** invokes Snakemake with `--forceall`, forcing all workflow rules to execute again using the current configuration. Existing output files may be overwritten if the rules complete successfully. The previous run's output files can be deleted fully by running:
```bash
rm -rf results/summary/* results/align/* results/raw/* results/qc/*
```
> **Warning:** This will permanently delete generated pipeline outputs. Do not run if you need to preserve previous results. Instead, save them to another folder before cleaning up.

---
## 🔜 Roadmap

**Note: Iostrack is a prototype and work in progress. Results should be independently reviewed before being used in publications, clinical workflows, or public-health decisions.
Features coming soon include: More customizability within plots, as well as more diverse and useful visualizations. More ways to query, filter, analyze and sort results. 
Any feedback is appreciated!** 

## 📝 License
This project is licensed under the MIT License. See [LICENSE](LICENSE).