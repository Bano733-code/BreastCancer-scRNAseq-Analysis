# 🧬 Single-Cell RNA-seq Analysis of the Breast Cancer Tumor Microenvironment

![R](https://img.shields.io/badge/R-4.6-276DC3?style=flat-square&logo=r)
![Seurat](https://img.shields.io/badge/Seurat-v5-blue?style=flat-square)
![CellChat](https://img.shields.io/badge/CellChat-v2-green?style=flat-square)
![License](https://img.shields.io/badge/License-MIT-lightgrey?style=flat-square)

An end-to-end scRNA-seq analysis of breast cancer (GSE176078) that goes beyond standard cell-type mapping: it classifies CNV-high (malignant-candidate) vs. CNV-low epithelial cells, characterizes malignant-state heterogeneity across clinical subtypes, infers cell-cell communication with CellChat (with a focus on midkine, MK), and evaluates a single-cell-derived gene signature as a prognostic marker in an independent bulk cohort (METABRIC).

<p align="center">
<img src="figures/Graphical_abstract/Graphical abstract.png" width="700">
</p>

---

## Dataset

| | |
|---|---|
| **Source** | GSE176078 — Wu et al., 2021, *Nat Genet* (26 primary, pre-treatment tumors) |
| **Cells** | 100,064 total → 94,195 after QC |
| **Samples** | 26 patients (11 ER+, 5 HER2+, 10 TNBC) |
| **Validation cohort** | METABRIC (n = 1,980), via cBioPortal |

The raw data are not part of this repository (see `.gitignore`): download GSE176078 into `data/GSE176078/GSE176078/` (including `metadata.csv`) and the METABRIC files into `data/metabric/metabric/`.

---

## Pipeline

```
QC → Normalization → PCA/UMAP → Clustering → SingleR annotation (+ comparison with the atlas labels)
        │
        ▼
Epithelial subset → chromosome-level CNV score → CNV-high / CNV-low
(cleaned-reference and held-out-null sensitivity analyses, benchmark against the atlas labels)
        │
        ▼
Malignant states → association with subtype (primary + all-patient sensitivity analysis)
        │
        ▼
CellChat (full object + subtype-stratified comparison) → MK checks (restricted network, patient level)
        │
        ▼
Gene signature → METABRIC scoring → survival analysis → signature benchmark (components, Ki-67, PAM50, random sets)
```

## Methods summary

- **QC & clustering:** Standard Seurat v5 workflow (cells with 300–7,000 detected genes and < 15% mitochondrial reads; 2,000 variable genes; 20 PCs; Louvain resolution 0.5); 31 clusters. Cell types annotated with SingleR (Human Primary Cell Atlas, `label.main`). About 10% of cells carry labels without a breast-tissue counterpart (e.g. tissue stem cells, chondrocytes); these are not interpreted biologically.
- **CNV-high classification:** A lightweight **chromosome-level** expression-deviation score (SD across chromosomes of the mean reference-adjusted expression; chromosome arms are *not* resolved; 23 chromosomes qualify). Immune and stromal cells form the reference; the cutoff is the 95th percentile of the reference scores (0.025). This proxy replaces inferCNV/CopyKAT for computational feasibility, so CNV-high cells are *malignant candidates* and CNV-low cells are *not confirmed normal*.
- **Malignant states:** Five curated states (basal-like, luminal progenitor-like, mature luminal-like, EMT/mesenchymal, cycling/proliferative) assigned with `AddModuleScore`; compared across subtypes with patient-level Kruskal-Wallis and Dunn tests.
- **Cell-cell communication:** CellChat (human CellChatDB, `triMean`, `raw.use = TRUE`) on the full object (default 100 permutations, ≤ 800 cells per group) and per subtype (25 permutations, ≤ 300 cells per group), with `liftCellChat` harmonization.
- **Signature & validation:** A 60-gene signature (50 cycling-state markers + 10 MK pathway genes; 55 present on the METABRIC array) scored as the mean gene-wise z-score; Kaplan-Meier and Cox models (OS, OS truncated at 120 months, disease-specific survival) with Schoenfeld tests.

---

## Key Results

| Metric | Value |
|---|---|
| CNV-high (malignant candidate) cells | 10,548 / 22,729 epithelial cells (46.4%) |
| Agreement with the atlas cancer/normal epithelial labels | specificity 99.9%, sensitivity 54.8% (original reference) |
| After removing cancer-like cells from the reference (10% contamination) | 16,561 CNV-high cells; sensitivity 86.0%, specificity 99.5%; held-out null gives the same cutoff (0.0209 vs 0.0211) |
| Malignant states associated with subtype | 5 / 5 in the primary analysis (BH p = 0.022–0.048); 4 / 5 when all 19 patients are kept (luminal progenitor-like p = 0.056) |
| Diversity of malignant states (patient-level Shannon index) | TNBC > ER+ (Dunn adj. p = 0.003) |
| CellChat pathways detected | 111 |
| Midkine (MK) relative information flow, TNBC / HER2+ / ER+ | 0.52 / 0.27 / 0.21 (full network); 0.40 / 0.31 / 0.29 (plausibly annotated groups only); higher MK potential in TNBC at the patient level (n = 16, p = 0.003) |
| Gene signature size | 60 genes (55 on the METABRIC array) |
| Survival HR, median split, METABRIC | 1.36 unadjusted; 1.31 (95% CI 1.12–1.53), p = 8.9 × 10⁻⁴ adjusted for age and subtype; 1.21 (1.03–1.43) fully adjusted |
| What carries the prognostic signal | the proliferation genes; the MK genes alone are not prognostic and add no significant information to them; the signature performs like the PAM50 proliferation score |

---

## Manuscript figures

All final panels are in `figures/Manuscript_Figures/` (scripts that regenerate them are listed in the table below).

**Figure 1 – Study overview and construction of the 60-gene signature**

<p align="center"><img src="figures/Manuscript_Figures/Fig1.png" width="760"></p>

**Figure 2 – Cellular landscape of the GSE176078 dataset** (a) clusters, (b) compartments, (c) cell counts, (d) marker genes

<table>
<tr><td><img src="figures/Manuscript_Figures/Fig2a.png" width="380"></td><td><img src="figures/Manuscript_Figures/Fig2b.png" width="380"></td></tr>
<tr><td><img src="figures/Manuscript_Figures/Fig2c.png" width="380"></td><td><img src="figures/Manuscript_Figures/Fig2d.png" width="380"></td></tr>
</table>

**Figure 3 – Identification of CNV-high epithelial cells** (a) epithelial cells on the original UMAP, (b) epithelial cells per sample, (c) distribution of chromosome-level CNV scores

<table>
<tr><td><img src="figures/Manuscript_Figures/Fig3a.png" width="300"></td><td><img src="figures/Manuscript_Figures/Fig3b.png" width="300"></td><td><img src="figures/Manuscript_Figures/Fig3c.png" width="300"></td></tr>
</table>

**Figure 4 – Malignant transcriptional-state heterogeneity** (a) UMAP of the 10,548 CNV-high cells, (b) module-score heatmap, (c) state composition by subtype, (d) patient-level proportions

<table>
<tr><td><img src="figures/Manuscript_Figures/Fig4a.png" width="380"></td><td><img src="figures/Manuscript_Figures/Fig4b.png" width="380"></td></tr>
<tr><td><img src="figures/Manuscript_Figures/Fig4c.png" width="380"></td><td><img src="figures/Manuscript_Figures/Fig4d.png" width="380"></td></tr>
</table>

**Figure 5 – Cell-cell communication and subtype-associated midkine signaling** (a) interaction strength, (b) relative information flow, (c) MK ligand-receptor bubble plot, (d) MK networks

<table>
<tr><td><img src="figures/Manuscript_Figures/Fig5a.jpg" width="380"></td><td><img src="figures/Manuscript_Figures/Fig5b.png" width="380"></td></tr>
<tr><td><img src="figures/Manuscript_Figures/Fig5c.png" width="380"></td><td><img src="figures/Manuscript_Figures/Fig5d.png" width="380"></td></tr>
</table>

**Figure 6 – Prognostic evaluation of the 60-gene signature** (a) score distribution in METABRIC, (b) Kaplan-Meier overall survival, (c) Cox regression

<table>
<tr><td><img src="figures/Manuscript_Figures/Fig6a.png" width="300"></td><td><img src="figures/Manuscript_Figures/Fig6b.jpg" width="300"></td><td><img src="figures/Manuscript_Figures/Fig6c.png" width="300"></td></tr>
</table>

### Which script makes which panel

| Panel | Script |
|---|---|
| Fig. 1 | schematic; the signature panel is produced by `figures/figures_scripts/Figure6(part A,B,D).R` (part 6A) |
| Fig. 2a | `figures/figures_scripts/Figure2(Cellular_Landscape).R` |
| Fig. 2b / 2c / 2d | `scripts/16_…`, `scripts/14_…`, `scripts/15_…` |
| Fig. 3a / 3b | `scripts/17_Epithelial_Subset+Reclustering.R` |
| Fig. 3c | `figures/figures_scripts/Figure3(part B,C).R` (last section) |
| Fig. 4a / 4b | `figures/figures_scripts/Figure4(Malignant_State_Module-Score_Heatmap).R` |
| Fig. 4c / 4d | `scripts/20_Malignant_Subtypes.R` |
| Fig. 5a | `scripts/21_CellChat.R` (`cellchat_overall_heatmap.pdf`) |
| Fig. 5b | `scripts/22_Subtype_Stratified_CellChat.R` (`cellchat_subtype_pathway_rank_FULL.pdf`) |
| Fig. 5c / 5d | `figures/figures_scripts/Figure5(MK_Signalling_Network).R` |
| Fig. 6a / 6c | `figures/figures_scripts/Figure6(part A,B,D).R` (parts 6B and 6D) |
| Fig. 6b | `scripts/25_Survival_Analysis.R` (`survival_KM_signature.pdf`) |

<p align="center">
<img src="results/TME/malignant_state_by_subtype_stackedbar.png" width="440">
<img src="results/TME/survival_KM_signature.jpg" width="440">
</p>

---

## Analyses added during peer review

Existing scripts were extended in place (every change is tagged `# [REVISION]`); only two new scripts were needed.

| Analysis | Where |
|---|---|
| Chromosome-level wording; reference-cell list; reference composition vs. atlas labels; cleaned-reference sensitivity analysis; held-out null; benchmark against the atlas labels; patient-level CNV table (Table 2) | `scripts/18_Malignant_vs_Normal_epithelial_cells.R` (sections 10–13) |
| SingleR vs. atlas label concordance | `scripts/08_Cell_Annotation.R` (end) |
| Cell counts per label and compartment (Table 1) | `scripts/16_Tumor_Microenvironment_Analysis.R` (end) |
| All-patient sensitivity analysis of malignant states (Dunn tests, leave-one-out, Shannon index) | `scripts/20_Malignant_Subtypes.R` (section 10) |
| CellChat settings documented | `scripts/21_CellChat.R`, `scripts/22_Subtype_Stratified_CellChat.R` |
| MK checks without re-running CellChat (restricted network, patient-level proxy) | `scripts/22_Subtype_Stratified_CellChat.R` (section 9) |
| Signature table S11 (duplicate symbols, renamed genes) | `scripts/24_Bulk_Validation.R` (end) |
| Extended survival models (DSS, 120-month truncation, per-SD hazard ratios, Schoenfeld tests) | `scripts/25_Survival_Analysis.R` (section 7) |
| Signature benchmark: proliferation vs. MK components, Ki-67, PAM50 proliferation score, 1,000 random gene sets | `scripts/26_Signature_benchmark.R` (new) |
| Supplementary workbook S1–S14 | `scripts/27_Supplementary_tables.R` (new) |

---

## Limitations

- CNV-high classification uses a simplified **chromosome-level** proxy, not inferCNV/CopyKAT. About 10% of the reference cells were cancer-epithelial according to the atlas labels, which made the original cutoff conservative: the CNV-high set is a high-purity but incomplete subset of the malignant cells, and CNV-low cells are not confirmed normal. A cleaned-reference sensitivity analysis is provided, but the downstream analyses use the original 10,548 cells.
- About 10% of cells carry SingleR labels without a breast-tissue counterpart; some of them (iPS- and MSC-labelled cells) are almost entirely cancer-epithelial cells according to the atlas.
- Subtype comparisons rest on few patients (7 TNBC, 8 ER+ and only 4 HER2+ patients contributed CNV-high cells); the association of luminal progenitor-like cells with subtype is not robust when all patients are kept, and the choice of a single "strongest" state is not robust.
- CellChat subtype analyses were run on downsampled cells with reduced permutations (`nboot = 25`); the MK difference between subtypes shrinks when groups with implausible labels are excluded.
- The gene signature is dominated by proliferation genes; the MK genes add no significant prognostic information beyond them, and the signature performs like the PAM50 proliferation score. Proportional hazards do not hold, tumor grade and size were not available, and only one bulk cohort (METABRIC) was used; TCGA-BRCA replication is planned.

---

## Repository Structure

```
├── figures/
│   ├── Manuscript_Figures/   # final manuscript panels (Fig1 … Fig6)
│   ├── figures_scripts/      # scripts that regenerate the manuscript figures
│   └── …                     # other analysis plots (QC, UMAP, GO/KEGG, malignant subtyping, CellChat, signature)
├── results/                  # tables and outputs (.csv/.pdf) underlying each figure
├── scripts/                  # analysis scripts 01–27, run sequentially
├── main.R
└── README.md
```

---

## Running the Pipeline

```r
install.packages(c("Seurat", "tidyverse", "patchwork", "data.table", "survival", "survminer",
                   "FSA", "pROC", "writexl"))
BiocManager::install(c("SingleR", "celldex", "clusterProfiler", "ReactomePA"))
devtools::install_github("sqjin/CellChat")
devtools::install_github("immunogenomics/presto")
```

Run the scripts in `scripts/` in numerical order from the project root (no `setwd()` calls are used). Random seed 123 is set wherever random sampling is used. The CNV scoring, CellChat and signature-benchmark steps are memory/time-intensive (see in-script notes for downsampling parameters).

---

## Supplementary data

The supplementary workbook (`Supplementary_Tables_S1-S14.xlsx`, built by `scripts/27_Supplementary_tables.R`) and the list of reference cells (`results/TME/reference_cells_script18.txt`) are deposited on Zenodo. Zenodo DOI: _to be added upon deposition_.

---

## Citation

- Wu, S.Z. et al. *A single-cell and spatially resolved atlas of human breast cancers.* Nat Genet 53, 1334–1347 (2021).
- Curtis, C. et al. *Nature* 486, 346–352 (2012) — METABRIC.
- Hao, Y. et al. *Nat Biotechnol* 42, 293–304 (2024) — Seurat v5.
- Jin, S. et al. *Nat Commun* 12, 1088 (2021) — CellChat.

---

## Author

**Bano Rani** — BS Bioinformatics, Department of Computer Science, University of Agriculture Faisalabad

Supervisor: Dr. Sumaira Nishat

⭐ If this project is useful, consider starring the repository.
