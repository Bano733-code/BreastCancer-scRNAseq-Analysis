# 🧬 Single-Cell RNA-seq Analysis of Breast Cancer Tumor Microenvironment: From Cellular Architecture to a Validated Prognostic Signature

![R](https://img.shields.io/badge/R-4.6-blue) ![Seurat](https://img.shields.io/badge/Seurat-5.x-orange) ![CellChat](https://img.shields.io/badge/CellChat-2.x-green) ![License](https://img.shields.io/badge/License-MIT-lightgrey) ![Status](https://img.shields.io/badge/Status-Active-brightgreen)

A comprehensive single-cell RNA sequencing analysis pipeline for breast cancer, extending beyond standard cell-type mapping to computationally distinguish malignant from normal epithelial cells, characterize malignant cell-state heterogeneity across molecular subtypes, infer cell-cell communication networks, and validate a derived gene signature as an independent prognostic biomarker in bulk clinical cohorts.

---

## Overview

This project performs an end-to-end scRNA-seq analysis of the breast cancer tumor microenvironment (TME), organized into four stages:

1. **Foundational pipeline** — QC, normalization, clustering, cell-type annotation, differential expression, and functional enrichment (31 clusters, immune/stromal/epithelial compartments).
2. **Phase 1 — Malignant Cell Subtyping** — CNV-based classification of malignant vs. normal epithelial cells, followed by malignant-state characterization and molecular-subtype association testing.
3. **Phase 2 — Cell-Cell Communication** — CellChat-based ligand-receptor signaling inference across the full annotated object, with subtype-stratified comparison.
4. **Phase 3 — Bulk Validation** — Derivation of a gene signature and its prognostic validation in an independent bulk RNA-seq cohort (METABRIC) via Kaplan-Meier and Cox regression survival analysis.

This design moves the project from *descriptive cellular mapping* (identifying that immune/stromal/epithelial populations exist) to a *mechanistic and translational result* (subtype-specific malignant heterogeneity, a candidate signaling axis, and a validated prognostic signature).

---

## Dataset

- **GEO Accession:** GSE176078
- **Organism:** Homo sapiens
- **Disease:** Breast Cancer
- **Technology:** Single-cell RNA Sequencing
- **Source:** Wu et al., 2021, *Nature Genetics* — "A single-cell and spatially resolved atlas of human breast cancers"
- **Scale:** 100,064 cells (26 patient samples); 94,195 retained after QC
- **Clinical metadata:** author-provided molecular subtype per sample (TNBC, HER2+, ER+), used throughout Phase 1 and Phase 2 for subtype-stratified analyses

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
┌─────────────────────────────────────────────┐
│  PHASE 1 — Malignant Cell Subtyping          │
│  Epithelial subset → CNV-based malignant/    │
│  normal classification → malignant states →  │
│  molecular subtype association                │
└─────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────┐
│  PHASE 2 — Cell-Cell Communication            │
│  CellChat (full object) → dominant pathway    │
│  identification → subtype-stratified          │
│  comparison                                    │
└─────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────┐
│  PHASE 3 — Bulk Validation                    │
│  Gene signature derivation → METABRIC         │
│  scoring → Kaplan-Meier + Cox survival        │
└─────────────────────────────────────────────┘
```

---

## Methods

### 1. Quality Control and Preprocessing
Cells filtered on: `nFeature_RNA > 300`, `nFeature_RNA < 7000`, mitochondrial percentage `< 15%`. 94,195 of 100,064 cells retained.

### 2. Dimensionality Reduction and Clustering
Highly variable genes → PCA → UMAP → unsupervised (Louvain) clustering, resolved into 31 clusters.

### 3. Cell Type Annotation
SingleR reference-based annotation identified T cells, NK cells, B cells, macrophages, monocytes, endothelial cells, fibroblasts, and epithelial cells, among others.

### 4. Differential Expression & Functional Enrichment
Cluster-specific marker genes (21,456 total DEGs; adj. p < 0.05, positive log2FC) were identified. GO Biological Process and Reactome pathway enrichment highlighted immune activation, antigen presentation, cytokine signaling, extracellular matrix organization, angiogenesis, and cancer-associated signaling.

### 5. Phase 1 — Malignant Cell Subtyping

**5.1 Epithelial subset & CNV-based classification.** Epithelial cells (SingleR call, n = 22,729) were extracted. Due to computational resource constraints, a lightweight chromosome-arm expression-deviation proxy (conceptually related to inferCNV) was implemented in place of inferCNV/CopyKAT: expression was referenced against a stratified sample of 4,743 immune/stromal cells, per-cell CNV scores were computed as the standard deviation of relative expression across chromosome arms, and a data-driven threshold (95th percentile of the reference distribution) classified cells as **CNV-high (malignant candidate)** or **CNV-low**. Result: **10,548 malignant candidate cells**, **12,181 CNV-low (normal) epithelial cells**; per-sample malignant fraction ranged 0–87%, consistent with expected inter-tumor heterogeneity.

**5.2 Malignant state characterization.** CNV-high cells were re-clustered independently (30 PCs, resolution 1.0) and annotated into five states via `AddModuleScore` using literature-curated marker sets: **Basal-like** (n=1,093), **Luminal Progenitor-like** (n=1,855), **Mature Luminal-like** (n=4,758), **EMT/Mesenchymal** (n=1,020), **Cycling/Proliferative** (n=1,822).

**5.3 Molecular subtype association.** Malignant state proportions were compared across TNBC/HER2+/ER+ at the patient (sample) level using Kruskal-Wallis tests with Benjamini-Hochberg correction, avoiding pseudoreplication from pooling correlated cells. **All five states showed significant subtype association** (adjusted p < 0.05), robust to exclusion of the largest-contributing sample. ER+ tumors were dominated by Mature Luminal-like cells (92.2%); HER2+ by Luminal Progenitor-like (50.5%); TNBC showed the greatest heterogeneity across all five states.

### 6. Phase 2 — Cell-Cell Communication (CellChat)

The full object was re-annotated with malignant-state labels replacing the generic "Epithelial" call. CellChat (human CellChatDB) was run on the full annotated object (111 active pathways detected, including the literature-expected TGF-β, SPP1, CXCL/CCL, and MIF axes) and separately on each molecular subtype (with cell-group harmonization via `liftCellChat` prior to merging/comparison). The unrestricted, full-database subtype comparison identified **midkine (MK) signaling** as the top differentially-active pathway between TNBC and HER2+/ER+ tumors — a data-driven finding rather than an assumed literature candidate.

### 7. Phase 3 — Bulk Validation & Survival Analysis

A 60-gene signature was derived by combining (a) the top 50 marker genes of the malignant state most strongly associated with subtype (Cycling/Proliferative), and (b) the 10 ligand-receptor genes of the MK signaling pathway. The signature was scored (mean z-score method) in the independent METABRIC bulk cohort (n = 1,980 patients; 55/60 genes present on the microarray platform) and patients stratified by median split. Kaplan-Meier and Cox proportional hazards regression (unadjusted and adjusted for age and molecular subtype) tested prognostic relevance.

**Result:** High-signature patients showed significantly worse overall survival (log-rank p = 1 × 10⁻⁷; adjusted Cox HR = 1.31, 95% CI 1.12–1.53, p = 8.9 × 10⁻⁴), independent of age and subtype.

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
├── data/
│   ├── GSE176078/                     # raw scRNA-seq data + author metadata
│   └── METABRIC/                      # bulk validation cohort (cBioPortal)
│
├── scripts/
│   ├── 01_Load_Data.R
│   ├── 02_QC.R
│   ├── 03_Filtering.R
│   ├── 04_Normalization.R
│   ├── 05_PCA.R
│   ├── 06_Clustering_UMAP.R
│   ├── 07_Marker_Genes.R
│   ├── 08_Cell_Annotation.R
│   ├── 09_DEG_Analysis.R
│   ├── 10_GO_KEGG.R
│   ├── 11_Reactome_Pathway.R
│   ├── 12_Marker_Heatmap.R
│   ├── 13_Marker_DotPlot.R
│   ├── 14_Cell_Composition.R
│   ├── 15_TME_Analysis.R
│   │
│   ├── phase1_malignant_subtyping/
│   │   ├── 16_CNV_Classification.R          # epithelial/reference selection + CNV scoring
│   │   ├── 17_Malignant_States.R             # re-clustering + AddModuleScore annotation
│   │   ├── 18_Subtype_Comparison.R           # chi-square/Fisher + stacked bar plots
│   │   └── 19_Subtype_Stats_Addendum.R       # sample-level Kruskal-Wallis + Dunn post-hoc
│   │
│   ├── phase2_cellchat/
│   │   ├── 20_Build_Annotated_Object.R       # merge malignant-state labels into full object
│   │   ├── 21_RunCellChat_Full.R             # full-object CellChat + pathway ID
│   │   └── 22_Subtype_Stratified_CellChat.R  # per-subtype CellChat + liftCellChat + rankNet
│   │
│   └── phase3_bulk_validation/
│       ├── 23_Gene_Signature.R               # signature derivation (state + pathway genes)
│       ├── 24_Bulk_Validation_METABRIC.R     # signature scoring in bulk cohort
│       └── 25_Survival_Analysis.R            # Kaplan-Meier + Cox regression
│
├── figures/
├── results/
│   └── TME/                                  # all intermediate .rds / .csv outputs
├── README.md
└── .gitignore
```

---

## Running the Pipeline

Install required packages:

```r
install.packages(c("Seurat", "tidyverse", "patchwork", "data.table",
                    "R.utils", "survival", "survminer"))

BiocManager::install(c("SingleR", "celldex", "clusterProfiler", "ReactomePA"))

devtools::install_github("sqjin/CellChat")
devtools::install_github("immunogenomics/presto")   # faster Wilcoxon test for CellChat
```

Run the foundational pipeline, then each phase in order:

```r
# Foundational pipeline
source("scripts/01_Load_Data.R")
...
source("scripts/15_TME_Analysis.R")

# Phase 1 — Malignant Cell Subtyping
source("scripts/phase1_malignant_subtyping/16_CNV_Classification.R")
source("scripts/phase1_malignant_subtyping/17_Malignant_States.R")
source("scripts/phase1_malignant_subtyping/18_Subtype_Comparison.R")
source("scripts/phase1_malignant_subtyping/19_Subtype_Stats_Addendum.R")

# Phase 2 — Cell-Cell Communication
source("scripts/phase2_cellchat/20_Build_Annotated_Object.R")
source("scripts/phase2_cellchat/21_RunCellChat_Full.R")
source("scripts/phase2_cellchat/22_Subtype_Stratified_CellChat.R")

# Phase 3 — Bulk Validation
source("scripts/phase3_bulk_validation/23_Gene_Signature.R")
source("scripts/phase3_bulk_validation/24_Bulk_Validation_METABRIC.R")
source("scripts/phase3_bulk_validation/25_Survival_Analysis.R")
```

**Note on compute resources:** Phase 1 CNV scoring and Phase 2 CellChat runs are memory- and time-intensive (CellChat's `computeCommunProb` can take 30–90+ minutes per run depending on cell/group counts). Scripts include downsampling and caching logic to keep runs feasible on a standard laptop (16 GB RAM). See in-script comments for tunable parameters (`max_cells_per_group`, `nboot`).

---

## Biological Insights

- Immune, stromal, and epithelial populations were resolved with characteristic markers (T cells: CD3D/CD3E/CD2/CCL5; NK cells: NKG7/GNLY; B cells: CD79A/MS4A1; macrophages: LST1/TYROBP; endothelial cells: VWF/EMCN/PLVAP; fibroblasts: COL1A1/COL3A1).
- Malignant epithelial cells are **not a homogeneous population**: five transcriptionally distinct states exist, and their relative proportions are significantly shaped by molecular subtype — TNBC tumors are markedly more heterogeneous than ER+ or HER2+ tumors.
- Midkine (MK) signaling emerges as a candidate subtype-differential communication axis, alongside established breast-TME pathways (TGF-β, SPP1, CXCL/CCL, MIF).
- A signature combining the top subtype-associated malignant state with the MK pathway gene set is independently prognostic for overall survival in an external cohort of nearly 2,000 patients.

---

## Limitations

This project prioritizes transparency about methodological trade-offs made under computational constraints:

- **CNV classification** used a simplified chromosome-arm expression-deviation proxy rather than inferCNV/CopyKAT; cells are reported as "CNV-high (malignant candidates)" rather than definitively malignant, and this method has not yet been benchmarked against gold-standard tools on this dataset.
- **Sample-level (patient-as-unit) statistics** were used as the primary evidence for subtype associations to avoid pseudoreplication inherent in treating individual cells as independent observations.
- **CellChat runs** used downsampled cell counts per group and reduced permutation counts (`nboot = 25` vs. the default 100) for computational feasibility.
- The final gene signature is dominated by canonical proliferation markers, which are near-universally prognostic in breast cancer cohorts; this result should be interpreted partly as a **positive-control validation of the overall pipeline**, and its value beyond established proliferation indices (e.g., Ki-67, PAM50 proliferation score) remains to be tested.
- Prognostic validation was performed in a **single bulk cohort (METABRIC)**; independent replication in TCGA-BRCA is a planned next step.

---

## 📌 Future Improvements

- [ ] Validate the gene signature in TCGA-BRCA as an independent second cohort
- [ ] Benchmark signature specificity against single-gene (Ki-67) and PAM50 proliferation-score baselines
- [ ] Subtype-stratified survival analysis (test whether the signature remains prognostic within TNBC specifically)
- [ ] Benchmark the CNV-score proxy method against inferCNV/CopyKAT
- [ ] Monocle3 trajectory inference across malignant states
- [ ] RNA velocity analysis
- [ ] Spatial transcriptomics integration
- [ ] Multi-omics integration

*(Completed since the original project scope: ~~CellChat analysis~~ ✅, ~~Copy number variation inference~~ ✅)*

---

## 📜 Citation

If this repository contributes to your research, please consider citing:

- The original dataset: Wu, S.Z. et al. *A single-cell and spatially resolved atlas of human breast cancers.* Nat Genet 53, 1334–1347 (2021).
- The METABRIC cohort: Curtis, C. et al. *Nature* 486, 346–352 (2012); Pereira, B. et al. *Nat Commun* 7, 11479 (2016).
- Software: Seurat (Hao et al., 2021), SingleR (Aran et al., 2019), CellChat (Jin et al., 2021).

---

## 👩‍💻 Author

**Bano Rani**
BS Bioinformatics
Department of Computer Science, University of Agriculture Faisalabad

**Research Interests:** Single-cell Genomics · Cancer Bioinformatics · Machine Learning · Computational Biology · AI for Precision Medicine

---

## ⭐ Support

If you find this project useful, please consider starring the repository.
