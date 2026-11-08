# ============================================================
# 27_Supplementary_tables.R
# Assembles the supplementary workbook "Supplementary_Tables_S1-S14.xlsx"
# (deposited on Zenodo) from the output files of scripts 08-26.
# Run AFTER all other scripts. Sheets whose source file name depends on your
# own earlier outputs are looked up below; a missing file only triggers a warning.
# ============================================================

if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
rd <- function(f) {
  if (!file.exists(f)) { warning("Missing file: ", f); return(data.frame(note = paste("file not found:", f))) }
  read.csv(f, check.names = FALSE)
}

# Cox results are split into median-split (S12) and continuous-score (S13) tables;
# p-values below 0.0001 for the Schoenfeld tests are shown as "<0.0001"
res <- read.csv("results/TME/survival_models_full.csv")
for (cn in c("PH_term_p", "PH_global_p")) res[[cn]] <- ifelse(res[[cn]] < 0.0001, "<0.0001", res[[cn]])

readme <- data.frame(
  Table = paste0("S", 1:14),
  Content = c(
    "Differential-expression genes; SingleR vs original atlas labels; SingleR x compartment counts",
    "Marker genes", "Pathway enrichment (Reactome)",
    "CNV scores per cell, cleaned-reference sensitivity analysis, held-out null, benchmark against atlas labels",
    "Malignant-state scores per cell",
    "Association of malignant states with subtype (primary analysis)",
    "Per-sample state proportions and sensitivity analyses (all patients, Dunn tests, leave-one-out, Shannon index)",
    "CellChat candidate pathways", "CellChat significant pathways",
    "CellChat subtype pathway ranking; midkine (MK) full-network, restricted-network and patient-level analyses",
    "60-gene signature (component, presence on the METABRIC array)",
    "Cox models, median split", "Cox models, continuous score",
    "Signature components, incremental likelihood-ratio tests, random gene sets, PAM50 proliferation score"))

sheets <- list(
  README = readme,
  S1_DEGs                 = rd("results/Differential_Expression/All_DEGs.csv"),
  S1_SingleR_vs_atlas     = rd("results/Annotation/SingleR_vs_atlas_crosstab.csv"),
  S1_SingleR_vs_TME       = rd("results/TME/SingleR_vs_TME_crosstab.csv"),
  S2_markers              = rd("results/Marker_Genes/All_Markers.csv"),
  S3_Reactome             = rd("results/Reactome/Reactome_Pathways.csv"),
  S4_cnv_per_cell         = rd("results/TME/cnv_scores_per_cell.csv"),
  S4_cleanRef_cells       = rd("results/TME/CNV_sensitivity_cleanRef_cells.csv"),
  S4_reference_atlas      = rd("results/TME/CNV_reference_composition_atlas.csv"),
  S4_heldout_null         = rd("results/TME/CNV_heldout_null.csv"),
  S4_benchmark            = rd("results/TME/CNV_benchmark_summary.csv"),
  S4_benchmark_patient    = rd("results/TME/CNV_benchmark_per_patient.csv"),
  S4_table2_patients      = rd("results/TME/Table2_per_patient_CNV.csv"),
  S5_state_scores         = rd("results/TME/malignant_state_scores_per_cell.csv"),
  S6_state_KW_primary     = rd("results/TME/malignant_state_kruskal_BHadjusted.csv"),
  S6_state_Dunn_primary   = rd("results/TME/malignant_state_dunn_posthoc.csv"),
  S6_state_association    = rd("results/TME/malignant_state_subtype_association.csv"),
  S7_sample_props         = rd("results/TME/sample_level_state_proportions.csv"),
  S7_all_patients         = rd("results/TME/malignant_state_per_patient.csv"),
  S7_KW_zero_incl         = rd("results/TME/S7_KW_zero_inclusive.csv"),
  S7_Dunn                 = rd("results/TME/S7_Dunn_zero_inclusive.csv"),
  S7_leave_one_out        = rd("results/TME/S7_leave_one_out.csv"),
  S7_shannon              = rd("results/TME/malignant_state_shannon.csv"),
  S8_candidate_pathways   = rd("results/TME/cellchat_candidate_pathways_table.csv"),
  S9_significant_pathways = rd("results/TME/cellchat_significant_pathways_sorted.csv"),
  S10_pathway_ranking     = rd("results/TME/cellchat_subtype_pathway_rank_table.csv"),
  S10_MK_flow_full        = rd("results/TME/MK_flow_full_network.csv"),
  S10_MK_flow_restricted  = rd("results/TME/MK_flow_restricted.csv"),
  S10_MK_per_patient      = rd("results/TME/MK_per_patient_restricted.csv"),
  S10_MK_leave_one_out    = rd("results/TME/MK_leave_one_patient_out_restricted.csv"),
  S11_signature           = rd("results/TME/Supplementary_Table_S11.csv"),
  S12_cox_median_split    = subset(res, term == "group"),
  S13_cox_continuous      = subset(res, term != "group"),
  S14_components          = rd("results/TME/signature_components.csv"),
  S14_incremental_LRT     = rd("results/TME/signature_incremental_LRT.csv"),
  S14_random_sets         = rd("results/TME/signature_random_benchmark.csv"),
  S14_PAM50               = rd("results/TME/signature_vs_PAM50_prolif.csv")
)
writexl::write_xlsx(sheets, "results/Supplementary_Tables_S1-S14.xlsx")
cat("Saved: results/Supplementary_Tables_S1-S14.xlsx with", length(sheets), "sheets\n")
