# Atopic Dermatitis Spatial Transcriptomics (Visium) Analysis Pipeline

## 📌 Overview
This repository contains an R-based analysis pipeline for integrating and clustering 10x Genomics Visium spatial transcriptomics data. The pipeline is optimized for **Seurat v5** and is designed to process multiple samples (19 samples) efficiently, avoiding common memory allocation and data structure errors during integration.

The dataset focuses on Atopic Dermatitis, comparing Lesional (LS), Non-lesional (NL), and Healthy Control (HC) skin tissues.

## 🛠 Prerequisites
The analysis is conducted in R and requires the following packages:
* `Seurat` (v5.0 or higher)
* `dplyr`
* `cowplot`
* `RColorBrewer`

## 🚀 Pipeline Workflow

### Part 1: Data Loading & Preprocessing
* Loaded 20 raw Visium samples.
* Performed Quality Control (QC). 1 sample was automatically excluded due to low quality.
* Log-normalization applied to the remaining 19 samples.
* Output: `SampleListVisium_LogNorm.rds`

### Part 2: Integration (`Visium_Integration_Part2.R`)
* **Feature Selection:** Selected 2,000 highly variable features across all 19 samples to optimize memory.
* **Dimensional Reduction:** Scaled data and ran PCA on each individual sample.
* **Integration Method:** Used Reciprocal PCA (RPCA) via `FindIntegrationAnchors` and `IntegrateData` to accurately batch-correct and merge the spatial datasets.
* Output: `IntegratedVisium_LogNorm.rds`

### Part 3: Clustering & Visualization (`Visium_Integration_Part3.R`)
* **Metadata Annotation:** Grouped samples into disease states (`LS`, `NL`, `HC`) using `case_when` based on sample IDs (`orig.ident`).
* **Dimensional Reduction:** Utilized the explicitly mapped variable features from the integrated assay for `ScaleData`, `RunPCA`, and `RunUMAP`.
* **Clustering:** Graph-based clustering using `FindNeighbors` and `FindClusters` (resolution = 0.5).
* **Visualization:** 
  * Generated UMAP plots to visualize disease group distributions.
  * Generated Grid Spatial plots (`SpatialDimPlot`) mapped across all 19 tissue slices using custom high-contrast color palettes.
* Output: `Final_IntegratedVisium.rds`, `UMAP_Visium_Integrated.png`, `SpatialClusterPlots_Visium_v5_CLR.png`

## 💡 Key Technical Fixes Included
This pipeline includes robust error-handling for common Seurat v5 integration issues:
1. **Missing `dimnames` Error Fix:** Prevented layer fragmentation errors during PCA by explicitly inheriting `VariableFeatures` from the integrated assay rather than re-running `FindVariableFeatures`.
2. **Memory Leak Prevention:** Optimized object passing and utilized RPCA for large-scale integration to avoid `bad_alloc` memory errors.
3. **Automated Metadata Grouping:** Safe string matching (Regex) applied to automatically label disease conditions regardless of sample input order.

## 📁 Directory Structure
Ensure your local environment matches the following directory structure before running the scripts:

```text
atopic_dermatitis_visium/
│
├── GSM5907096/
│   ├── rds_OUT/                  # Saved .rds R objects
│   ├── visualization/
│   │   └── plots/                # Output directory for UMAP and Spatial PNGs
│   ├── Visium_Integration_Part2.R
│   └── Visium_Integration_Part3.R
