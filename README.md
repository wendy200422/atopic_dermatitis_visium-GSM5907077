# 🧬 Atopic Dermatitis Transcriptomics Analysis Pipeline

## 📌 Overview
This repository contains R-based analysis pipelines for **Spatial Transcriptomics (10x Visium)** and **Single-cell RNA-seq (scRNA-seq)** data. The workflow focuses on analyzing Atopic Dermatitis skin tissues, specifically comparing Lesional (LS), Non-lesional (NL), and Healthy Control (HC) states.

* **Visium Dataset:** Sourced from GEO under accession number **[GSE197023](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197023)**.
* **Key Features:** Optimized for **Seurat v5**, addressing memory allocation issues and ensuring stable dimensional reduction and clustering in dry-lab environments.

## 📚 References & Acknowledgements
* **Reference Pipeline:** The spatial transcriptomics workflow was closely adapted from and inspired by [juno-kim/visium-atopic-dermatitis](https://gitlab.com/juno-kim/visium-atopic-dermatitis).

## 📁 Repository Structure
```text
script/
├── 📄 README.md
├── 🔬 Visium (Spatial Transcriptomics)
│   ├── Visium_Integration_Part1_v5.R
│   ├── Visium_Integration_Part2_v5.R
│   └── Visium_Integration_Part3_v5.R
└── 🧫 scRNA-seq (Single-cell RNA-seq)
    ├── scRNA-seq_Part1.R
    ├── scRNA-seq_Part2.R
    ├── scRNA-seq_Part3.R
    ├── scRNA-seq_Part4.R
    ├── scRNA-seq_Part5.R
    └── scRNA-seq_DGE.R
