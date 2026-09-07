# 🧬 Single-Cell RNA-seq Analysis of Breast Cancer Tumor Microenvironment: From Cellular Architecture to a Validated Prognostic Signature

<div align="center">

![R](https://img.shields.io/badge/R-4.6-276DC3?style=for-the-badge&logo=r)
![Seurat](https://img.shields.io/badge/Seurat-v5-blue?style=for-the-badge)
![CellChat](https://img.shields.io/badge/CellChat-v2-green?style=for-the-badge)
![License](https://img.shields.io/badge/License-MIT-lightgrey?style=for-the-badge)
![Status](https://img.shields.io/badge/Project-Completed-success?style=for-the-badge)

*A comprehensive single-cell RNA sequencing analysis pipeline for breast cancer, extending beyond standard cell-type mapping to malignant-cell classification, subtype-specific heterogeneity, cell-cell communication inference, and a validated prognostic gene signature.*

</div>

---

## Overview

This project performs an end-to-end scRNA-seq analysis of the breast cancer tumor microenvironment (TME), organized into four stages:

1. **Foundational pipeline** — QC, normalization, clustering, cell-type annotation, differential expression, and functional enrichment.
2. **Malignant cell subtyping** — CNV-based classification of malignant candidate vs. CNV-low epithelial cells, followed by malignant-state characterization and clinical-subtype association testing.
3. **Cell-cell communication** — CellChat-based ligand-receptor signaling inference across the full annotated object, with subtype-stratified comparison.
4. **Bulk validation** — Derivation of a gene signature and its prognostic validation in an independent bulk RNA-seq cohort (METABRIC) via Kaplan-Meier and Cox regression survival analysis.

---

## Dataset

- **GEO Accession:** GSE176078
- **Organism:** Homo sapiens
- **Disease:** Breast Cancer
- **Technology:** Single-cell RNA Sequencing
- **Source:** Wu et al., 2021, *Nature Genetics* — "A single-cell and spatially resolved atlas of human breast cancers"
- **Scale:** 100,064 cells (26 patient samples); 94,195 retained after QC
- **Clinical metadata:** author-provided clinical subtype per sample (11 ER+, 5 HER2+, 10 TNBC)

---

## Workflow

```
Raw scRNA-seq Data
        │
        ▼
QC, Filtering, Normalization
        │
        ▼
PCA → UMAP → Clustering (31 clusters)
        │
        ▼
Cell Type Annotation (SingleR)
        │
        ▼
Differential Expression + GO / Reactome Enrichment
        │
        ▼
Epithelial Subset → CNV-Based Classification → Malignant States
        │                                            │
        │                                            ▼
        │                              Clinical Subtype Association
        ▼
Full Annotated Object → CellChat → Subtype-Stratified Comparison
        │
        ▼
Gene Signature Derivation → METABRIC Scoring → Survival Analysis
```

---

## 1. Quality Control and Preprocessing

Cells were filtered on:
```
nFeature_RNA > 300
nFeature_RNA < 7000
Mitochondrial percentage < 15%
```
94,195 of 100,064 cells were retained.

<p align="center">
<img src="figures/QC/QC_Violin.png" width="700"><br>
<img src="figures/QC/QC_Scatter.png" width="700"><br>
</p>

---

## 2. Dimensionality Reduction and Clustering

Highly variable genes were identified, followed by PCA and UMAP-based visualization of unsupervised (Louvain) clustering, resolving 31 clusters.

<p align="center">
<img src="figures/UMAP/UMAP.png" width="700">
</p>

---

## 3. Cell Type Annotation

Cell identities were assigned using **SingleR** reference-based annotation:

| Cell Population | Description |
|---|---|
| T Cells | Adaptive immune population |
| NK Cells | Cytotoxic lymphocytes |
| B Cells | Antibody-producing immune cells |
| Macrophages | Tumor-associated immune cells |
| Monocytes | Myeloid immune population |
| Endothelial Cells | Tumor vasculature |
| Fibroblasts | Stromal population |
| Epithelial Cells | Tumor-associated cells |

<p align="center">
<img src="figures/Marker_DotPlot/Marker_DotPlot.png" width="700">
</p>

**Representative markers:**

| Cell Type | Markers |
|---|---|
| T Cells | CD3D, CD3E, CD2, CCL5 |
| NK Cells | NKG7, GNLY |
| B Cells | CD79A, MS4A1 |
| Macrophages | LST1, TYROBP |
| Endothelial Cells | VWF, EMCN, PLVAP |
| Fibroblasts | COL1A1, COL3A1 |

---

## 4. Differential Expression & Functional Enrichment

Cluster-specific DEGs (21,456 total; adjusted p < 0.05, positive log2FC) were identified. GO Biological Process and Reactome pathway enrichment highlighted immune activation, antigen presentation, cytokine signaling, extracellular matrix organization, and angiogenesis.

<p align="center">
<img src="figures/GO_KEGG/GO_Dotplot.png" width="700"><br>
<img src="figures/GO_KEGG/GO_Barplot.png" width="700">
</p>

---

## 5. Tumor Microenvironment Composition

<p align="center">
<img src="figures/TME/TME_UMAP.png" width="700"><br>
<img src="figures/TME/Immune_Landscape.png" width="700">
</p>

| Compartment | Cell Types |
|---|---|
| Immune | T cells, NK cells, B cells, Macrophages, Monocytes |
| Stromal | Fibroblasts, Endothelial cells |
| Tumor/Epithelial | Epithelial cells (subsequently reclassified — see Section 6) |

<p align="center">
<img src="figures/Cell_Composition/Cell_Composition_Barplot.png" width="600">
<img src="figures/Cell_Composition/Cell_Composition_PieChart.png" width="600">
</p>

---

## 6. Malignant Cell Subtyping

### 6.1 Epithelial re-clustering
Epithelial cells (n = 22,729) were extracted and re-normalized/re-clustered independently.

<p align="center">
<img src="figures/Malignant/Epithelial_Cells_Original_UMAP.png" width="600">
<img src="figures/Malignant/Epithelial_Subclusters_UMAP.png" width="600"><br>
<img src="figures/Malignant/Epithelial_PCA_Elbow.png" width="500">
<img src="figures/Malignant/Epithelial_Variable_Features.png" width="500">
</p>

### 6.2 CNV-based malignant classification
Due to computational constraints, a lightweight chromosome-arm expression-deviation proxy (conceptually related to inferCNV) was used in place of inferCNV/CopyKAT. Expression was referenced against a stratified sample of 4,743 immune/stromal cells; a data-driven threshold (95th percentile of the reference distribution) classified epithelial cells as **CNV-high (malignant candidate)** or **CNV-low**.

**Result:** 10,548 malignant candidate cells (46.4%); 12,181 CNV-low cells. Per-sample malignant fraction ranged 0–87%.

<p align="center">
<img src="figures/Malignant/Figure3B_CNV_Score_Distribution.png" width="600">
<img src="figures/Malignant/Figure3C_CNV_High_vs_Low.png" width="600"><br>
<img src="figures/Malignant/Epithelial_Cells_By_Sample.png" width="700">
</p>

### 6.3 Malignant state characterization
CNV-high cells were re-clustered and annotated into five states via `AddModuleScore`:

| State | n cells | Representative markers |
|---|---|---|
| Basal-like | 1,093 | KRT5, KRT14, TP63 |
| Luminal Progenitor-like | 1,855 | KIT, ELF5, ALDH1A3 |
| Mature Luminal-like | 4,758 | ESR1, FOXA1, GATA3 |
| EMT/Mesenchymal | 1,020 | VIM, ZEB1, FN1 |
| Cycling/Proliferative | 1,822 | MKI67, TOP2A, CDK1 |

<p align="center">
<img src="figures/Malignant/Figure4B_Malignant_State_ModuleScore_Heatmap.png" width="700">
</p>

### 6.4 Clinical subtype association
Malignant state proportions were compared across TNBC/HER2+/ER+ using patient-level Kruskal-Wallis testing with Benjamini-Hochberg correction (primary evidence, to avoid pseudoreplication from pooling correlated cells).

**All five states showed significant subtype association** (adjusted p = 0.022–0.048). ER+ tumors were dominated by Mature Luminal-like cells (92.2%); HER2+ by Luminal Progenitor-like (50.5%); TNBC displayed the broadest distribution across all five states.

<p align="center">
<img src="results/TME/malignant_state_by_subtype_stackedbar.png" width="600">
<img src="results/TME/malignant_state_by_subtype_persample_boxplot.png" width="600">
</p>

---

## 7. Cell-Cell Communication (CellChat)

The full object was re-annotated with malignant-state labels replacing the generic epithelial call. CellChat (human CellChatDB) was run on the full annotated object and separately per clinical subtype (with `liftCellChat` harmonization prior to merging).

**111 active signaling pathways** were detected, including the literature-expected TGF-β, SPP1, CXCL/CCL, and MIF axes. An unrestricted, full-database subtype comparison identified **midkine (MK) signaling** as the top subtype-differential pathway between TNBC and HER2+/ER+ tumors — a data-driven finding, not an assumed candidate.

<p align="center">
<img src="figures/TME/Figure5D_MK_Signaling_Network.png" width="700">
</p>

*Additional CellChat outputs (bubble plots, pathway ranking, interaction-strength comparisons) are in `results/TME/` — see `cellchat_subtype_pathway_rank_table.csv` for the full ranked pathway list and `cellchat_candidate_pathways_table.csv` for the five literature-expected pathways across subtypes.*

---

## 8. Gene Signature & Bulk Validation

A 60-gene signature was derived by combining the top 50 marker genes of the malignant state most strongly associated with subtype (Cycling/Proliferative) with the 10 ligand-receptor genes of the MK pathway (MDK, SDC1, SDC2, SDC4, PTPRZ1, LRP1, NCL, ITGA4, ITGA6, ITGB1).

<p align="center">
<img src="figures/Signature/Figure6A_60Gene_Signature_Construction_page-0001.jpg" width="700">
</p>

The signature was scored (mean z-score method) in the independent **METABRIC** cohort (n = 1,980; 55/60 genes present on the array platform) and patients stratified by median split.

<p align="center">
<img src="figures/Signature/Figure6B_METABRIC_Signature_Score_Distribution.png" width="600">
</p>

---

## 9. Survival Analysis

Kaplan-Meier and Cox proportional hazards regression (unadjusted and adjusted for age and clinical subtype) tested prognostic relevance.

**Result:** High-signature patients showed significantly worse overall survival (log-rank p = 1 × 10⁻⁷; adjusted Cox HR = 1.31, 95% CI 1.12–1.53, p = 8.9 × 10⁻⁴), independent of age and subtype. The continuous score (HR = 1.26, p = 1.9 × 10⁻⁷) confirmed this was not an artifact of the median-split cutpoint.

<p align="center">
<img src="figures/Signature/Figure6D_Publication_Cox_Plot_page-0001.jpg" width="700">
</p>

*Full survival outputs (KM curve, Cox summary table) are in `results/TME/survival_KM_signature.pdf` and `results/TME/survival_cox_summary.csv`.*

---

## Key Results Summary

| Metric | Value |
|---|---|
| Total cells (QC-passed) | 94,195 |
| Clusters | 31 |
| Epithelial cells | 22,729 |
| Malignant candidate cells (CNV-high) | 10,548 (46.4%) |
| Malignant states identified | 5 |
| States with significant subtype association | 5 / 5 (BH-adj. p < 0.05) |
| CellChat pathways detected | 111 |
| Top subtype-differential pathway | MK (midkine) |
| Final gene signature size | 60 genes |
| Validation cohort | METABRIC (n = 1,980) |
| Survival HR (adjusted) | 1.31 (95% CI 1.12–1.53), p = 8.9 × 10⁻⁴ |

---

## Repository Structure

```
BreastCancer-scRNAseq-Analysis/
│
├── figures/
├── results/
│   ├── Cell_Composition/
│   ├── Differential_Expression/
│   ├── GO_KEGG/
│   ├── Immune_Landscape/
│   ├── Malignant/
│   ├── Malignant_Subtyping/
│   ├── Marker_Genes/
│   ├── Reactome/
│   └── TME/  
│
├── scripts/
├── .gitignore
├── README.md
└── main.R
```

---

## Running the Pipeline

```r
install.packages(c("Seurat", "tidyverse", "patchwork", "data.table",
                    "R.utils", "survival", "survminer"))

BiocManager::install(c("SingleR", "celldex", "clusterProfiler", "ReactomePA"))

devtools::install_github("sqjin/CellChat")
devtools::install_github("immunogenomics/presto")   # faster Wilcoxon test for CellChat
```

Scripts are run sequentially from `scripts/` (foundational pipeline: QC → clustering → annotation → DE/enrichment), followed by the malignant subtyping, CellChat, and bulk validation stages described above.

**Note on compute resources:** CNV scoring and CellChat runs are memory- and time-intensive (CellChat's `computeCommunProb` can take 30–90+ minutes per run). Scripts include downsampling and caching logic to keep runs feasible on a standard laptop (16 GB RAM).

---

## Limitations

- **CNV classification** used a simplified chromosome-arm expression-deviation proxy rather than inferCNV/CopyKAT; cells are reported as "CNV-high (malignant candidates)" rather than definitively malignant, and this method has not yet been benchmarked against gold-standard tools on this dataset.
- **Sample-level (patient-as-unit) statistics** were used as the primary evidence for subtype associations to avoid pseudoreplication.
- **CellChat runs** used downsampled cell counts and reduced permutation counts (`nboot = 25` vs. the default 100) for computational feasibility.
- The gene signature is dominated by canonical proliferation markers, near-universally prognostic in breast cancer cohorts; this result is best read as a **positive-control validation of the pipeline**, and its value beyond established proliferation indices (e.g., Ki-67, PAM50) remains to be tested.
- Prognostic validation was performed in a **single bulk cohort (METABRIC)**; independent replication in TCGA-BRCA is a planned next step.

## 📌 Future Improvements

- [ ] Validate the gene signature in TCGA-BRCA as an independent cohort
- [ ] Decompose the signature into proliferation-only vs. MK-only components to test whether MK genes add prognostic value
- [ ] Benchmark signature specificity against Ki-67 / PAM50 proliferation-score baselines
- [ ] Subtype-stratified survival analysis (test whether the signature remains prognostic within TNBC specifically)
- [ ] Benchmark the CNV-score proxy method against inferCNV/CopyKAT
- [ ] Monocle3 trajectory inference across malignant states
- [ ] RNA velocity analysis
- [ ] Spatial transcriptomics integration

*(Completed since the original project scope: ~~CellChat analysis~~ ✅ ~~Copy number variation inference~~ ✅)*

---

## 📜 Citation

- Dataset: Wu, S.Z. et al. *A single-cell and spatially resolved atlas of human breast cancers.* Nat Genet 53, 1334–1347 (2021).
- METABRIC: Curtis, C. et al. *Nature* 486, 346–352 (2012); Pereira, B. et al. *Nat Commun* 7, 11479 (2016).
- Software: Seurat v5 (Hao et al., *Nat Biotechnol* 2024); SingleR (Aran et al., *Nat Immunol* 2019); CellChat (Jin et al., *Nat Commun* 2021).

---

## 👩‍💻 Author

**Bano Rani**
BS Bioinformatics · Department of Computer Science, University of Agriculture Faisalabad

**Supervisor:** Dr. Sumaira Nishat, Department of Computer Science, University of Agriculture Faisalabad

**Research Interests:** Single-cell Genomics · Cancer Bioinformatics · Machine Learning · Computational Biology · AI for Precision Medicine

---

## ⭐ Support

If you find this project useful, please consider starring the repository.
