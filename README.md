<p align="center">
<img src="figures/banner.png" width="100%">
</p>

# 🧬 Breast Cancer Single-Cell RNA-seq Analysis using Seurat

<div align="center">

![R](https://img.shields.io/badge/R-4.6+-276DC3?style=for-the-badge&logo=r)
![Seurat](https://img.shields.io/badge/Seurat-v5-blue?style=for-the-badge)
![SingleR](https://img.shields.io/badge/SingleR-Annotation-green?style=for-the-badge)
![License](https://img.shields.io/badge/License-MIT-red?style=for-the-badge)
![Status](https://img.shields.io/badge/Project-Completed-success?style=for-the-badge)

*A comprehensive single-cell RNA sequencing analysis pipeline for breast cancer using Seurat, SingleR, Gene Ontology, KEGG, Reactome, and Tumor Microenvironment characterization.*

</div>

---

# 📖 Overview

This repository presents a complete **single-cell RNA sequencing (scRNA-seq)** analysis workflow for breast cancer using the **GSE176078** dataset.

The project performs end-to-end analysis beginning with raw gene expression matrices and ending with biological interpretation through:

- Quality Control
- Normalization
- Dimensionality Reduction
- Clustering
- Cell Type Annotation
- Differential Expression Analysis
- Functional Enrichment
- Tumor Microenvironment Analysis
- Immune Landscape Characterization

The workflow follows current best practices commonly used in computational biology and cancer genomics research.

---

# 🎯 Objectives

The primary objectives of this project are:

- Analyze breast cancer single-cell transcriptomic data.
- Identify distinct cellular populations.
- Annotate cell types using automated reference-based methods.
- Discover cluster-specific marker genes.
- Characterize the tumor microenvironment.
- Investigate biological pathways enriched in different cell populations.
- Produce publication-quality visualizations.

---

# 🧬 Dataset

**Dataset**

GSE176078

**Species**

Homo sapiens

**Disease**

Breast Cancer

**Platform**

Single-cell RNA Sequencing

**Input Files**

- matrix.mtx
- genes.tsv
- barcodes.tsv
- metadata.csv

---

# 🔬 Analysis Workflow

```
Raw Expression Matrix
        │
        ▼
Quality Control
        │
        ▼
Filtering
        │
        ▼
Normalization
        │
        ▼
Highly Variable Genes
        │
        ▼
Scaling
        │
        ▼
Principal Component Analysis
        │
        ▼
Elbow Plot
        │
        ▼
UMAP
        │
        ▼
Clustering
        │
        ▼
Marker Gene Detection
        │
        ▼
SingleR Cell Annotation
        │
        ▼
Differential Expression
        │
        ▼
GO Enrichment
        │
        ▼
KEGG Pathway Analysis
        │
        ▼
Reactome Analysis
        │
        ▼
Tumor Microenvironment Analysis
        │
        ▼
Immune Landscape
```

---

# 📂 Project Structure

```
BreastCancer-scRNAseq/

│
├── data/
│   ├── raw/
│   ├── processed/
│
├── scripts/
│   ├── 01_Load_Data.R
│   ├── 02_Quality_Control.R
│   ├── 03_Filtering.R
│   ├── 04_Normalization.R
│   ├── 05_PCA.R
│   ├── 06_Clustering_UMAP.R
│   ├── 07_Marker_Genes.R
│   ├── 08_Cell_Annotation.R
│   ├── 09_Differential_Expression.R
│   ├── 10_GO_KEGG.R
│   ├── 11_Reactome_Pathway.R
│   ├── 12_Marker_Heatmap.R
│   ├── 13_Marker_DotPlot.R
│   ├── 14_Cell_Composition_and_Immune_Landscape.R
│   ├── 15_Tumor_Microenvironment_Analysis.R
│
├── figures/
│
├── results/
│
├── README.md
├── requirements.txt
├── LICENSE
└── .gitignore
```

---

# 📊 Analysis Performed

| Analysis | Status |
|-----------|---------|
| Data Loading | ✅ |
| Quality Control | ✅ |
| Cell Filtering | ✅ |
| Normalization | ✅ |
| Variable Gene Selection | ✅ |
| Scaling | ✅ |
| PCA | ✅ |
| UMAP | ✅ |
| Clustering | ✅ |
| Marker Gene Detection | ✅ |
| Cell Annotation (SingleR) | ✅ |
| Differential Expression | ✅ |
| GO Enrichment | ✅ |
| KEGG Analysis | ✅ |
| Reactome Pathway Analysis | ✅ |
| Marker Heatmap | ✅ |
| Marker DotPlot | ✅ |
| Cell Composition | ✅ |
| Immune Landscape | ✅ |
| Tumor Microenvironment Analysis | ✅ |

---

# 📈 Results

This analysis successfully identified multiple biologically meaningful cell populations, including:

- T Cells
- B Cells
- NK Cells
- Macrophages
- Monocytes
- Endothelial Cells
- Fibroblasts
- Epithelial (Tumor) Cells

Cluster-specific marker genes were identified and subsequently analyzed using GO, KEGG, and Reactome pathway enrichment to reveal functional characteristics of each cell population.

---

# 🖼 Example Figures

## Quality Control

<p align="center">
<img src="figures/QC/QC_Violin.png" width="700">
</p>

---

## PCA

<p align="center">
<img src="figures/PCA/PCA.png" width="700">
</p>

---

## UMAP

<p align="center">
<img src="figures/UMAP/UMAP.png" width="700">
</p>

---

## Marker Heatmap

<p align="center">
<img src="figures/Markers/Marker_Heatmap.png" width="700">
</p>

---

## GO Enrichment

<p align="center">
<img src="figures/GO/GO_BP.png" width="700">
</p>

---

## KEGG Pathway Analysis

<p align="center">
<img src="figures/KEGG/KEGG.png" width="700">
</p>

---

## Tumor Microenvironment

<p align="center">
<img src="figures/TME/TME_UMAP.png" width="700">
</p>

---

# 💻 Technologies Used

- R
- Seurat
- SingleR
- celldex
- clusterProfiler
- enrichplot
- ReactomePA
- ggplot2
- dplyr
- patchwork

---

# 🚀 Installation

```r
install.packages("Seurat")
install.packages("patchwork")
install.packages("tidyverse")

BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("clusterProfiler")
BiocManager::install("ReactomePA")
```

---

# ▶️ Running the Pipeline

Execute scripts sequentially.

```
01_Load_Data.R

↓

02_Quality_Control.R

↓

03_Filtering.R

↓

04_Normalization.R

↓

05_PCA.R

↓

06_Clustering_UMAP.R

↓

07_Marker_Genes.R

↓

08_Cell_Annotation.R

↓

09_Differential_Expression.R

↓

10_GO_KEGG.R

↓

11_Reactome_Pathway.R

↓

12_Marker_Heatmap.R

↓

13_Marker_DotPlot.R

↓

14_Cell_Composition_and_Immune_Landscape.R

↓

15_Tumor_Microenvironment_Analysis.R
```

---

# 📚 Biological Insights

Key findings include:

- Identification of multiple immune and stromal populations.
- Characterization of the breast cancer tumor microenvironment.
- Cluster-specific marker genes associated with immune activation.
- Functional enrichment of pathways related to immune response, extracellular matrix organization, angiogenesis, and cell proliferation.

---

# 📌 Future Improvements

Potential future extensions include:

- CellChat analysis for cell-cell communication
- Monocle3 trajectory inference
- RNA velocity analysis
- Copy number variation inference
- Spatial transcriptomics integration
- Multi-omics integration

---

# 📜 Citation

If this repository contributes to your research, please consider citing the original GSE176078 dataset and the Seurat framework.

---

# 👩‍💻 Author

**Bano Rani**

BS Bioinformatics

University of Agriculture Faisalabad

Research Interests:

- Single-cell Genomics
- Cancer Bioinformatics
- Machine Learning
- Computational Biology
- AI for Precision Medicine

---

# ⭐ Support

If you find this project useful, consider giving it a ⭐ on GitHub.
