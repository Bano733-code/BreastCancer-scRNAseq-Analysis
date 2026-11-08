# ============================================================
# Subtype-Stratified CellChat Comparison (TNBC / HER2+ / ER+)
# ============================================================

library(Seurat)
library(CellChat)
library(patchwork)
library(dplyr)

set.seed(123)

candidate_pathways <- c("TGFb", "SPP1", "CXCL", "CCL", "MIF")


# ------------------------------------------------------------
# 1. Reusable function: run CellChat for one subtype
#    (only called if the cached .rds doesn't already exist)
# ------------------------------------------------------------

run_cellchat_subtype <- function(seurat_obj, subtype_label,
                                 min_group_size = 15,
                                 max_cells_per_group = 300,
                                 nboot = 25) {
  
  cat("\nRunning CellChat for subtype:", subtype_label, "\n")
  
  sub_obj <- subset(seurat_obj, subset = subtype == subtype_label)
  
  celltype_counts <- table(sub_obj$final_celltype)
  keep_types <- names(celltype_counts)[celltype_counts >= min_group_size]
  sub_obj <- subset(sub_obj, subset = final_celltype %in% keep_types)
  
  set.seed(123)
  cells_keep <- unlist(lapply(unique(sub_obj$final_celltype), function(ct) {
    cells_in_type <- colnames(sub_obj)[sub_obj$final_celltype == ct]
    if (length(cells_in_type) > max_cells_per_group) {
      sample(cells_in_type, max_cells_per_group)
    } else {
      cells_in_type
    }
  }))
  sub_obj <- subset(sub_obj, cells = cells_keep)
  
  data_input <- GetAssayData(sub_obj, assay = "RNA", layer = "data")
  meta <- data.frame(labels = sub_obj$final_celltype, row.names = colnames(sub_obj))
  rm(sub_obj); gc()
  
  cc <- createCellChat(object = data_input, meta = meta, group.by = "labels")
  cc@DB <- CellChatDB.human
  
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  # [REVISION] Subtype-stratified analyses use 25 permutations (default is 100),
  # groups < 15 cells are dropped and groups are downsampled to <= 300 cells
  # to keep the run memory-safe; default 'triMean' method and raw.use = TRUE.
  cc <- computeCommunProb(cc, raw.use = TRUE, nboot = nboot)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  
  gc()
  return(cc)
}


# ------------------------------------------------------------
# 2. Load cached results if present, otherwise compute fresh
# ------------------------------------------------------------

get_or_compute <- function(path, subtype_label, seurat_obj) {
  if (file.exists(path)) {
    cat("Loading cached:", path, "\n")
    return(readRDS(path))
  } else {
    cc <- run_cellchat_subtype(seurat_obj, subtype_label)
    saveRDS(cc, path)
    return(cc)
  }
}

paths <- c(
  TNBC = "results/TME/cellchat_TNBC.rds",
  HER2 = "results/TME/cellchat_HER2.rds",
  ER   = "results/TME/cellchat_ER.rds"
)

if (!all(file.exists(paths))) {
  seurat_obj <- readRDS("results/TME/full_annotated_for_cellchat.rds")
  subtype_meta <- read.csv("data/GSE176078/GSE176078/metadata.csv")
  sample_subtype <- unique(subtype_meta[, c("orig.ident", "subtype")])
  seurat_obj$subtype <- sample_subtype$subtype[
    match(seurat_obj$orig.ident, sample_subtype$orig.ident)
  ]
} else {
  seurat_obj <- NULL   # not needed — everything cached
}

cellchat_TNBC <- get_or_compute(paths["TNBC"], "TNBC", seurat_obj)
cellchat_HER2 <- get_or_compute(paths["HER2"], "HER2+", seurat_obj)
cellchat_ER   <- get_or_compute(paths["ER"],   "ER+",   seurat_obj)

rm(seurat_obj); gc()


# ------------------------------------------------------------
# 3. Harmonize cell-type groups across all 3 objects (THE FIX)
#    — required before merging/comparing, since the three
#    subtypes don't contain identical sets of rare cell types
# ------------------------------------------------------------

group.new <- Reduce(union, list(
  levels(cellchat_TNBC@idents),
  levels(cellchat_HER2@idents),
  levels(cellchat_ER@idents)
))

cellchat_TNBC <- liftCellChat(cellchat_TNBC, group.new)
cellchat_HER2 <- liftCellChat(cellchat_HER2, group.new)
cellchat_ER   <- liftCellChat(cellchat_ER, group.new)


# ------------------------------------------------------------
# 4. Merge
# ------------------------------------------------------------

object.list <- list(TNBC = cellchat_TNBC, HER2 = cellchat_HER2, ER = cellchat_ER)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

saveRDS(cellchat_merged, "results/TME/cellchat_merged_comparison_FIXED.rds")


# ------------------------------------------------------------
# 5. Overall interaction strength (bar chart)
# ------------------------------------------------------------

pdf("results/TME/cellchat_subtype_interaction_strength.pdf", width = 8, height = 5)
print(compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2, 3)))
dev.off()


# ------------------------------------------------------------
# 6. rankNet — full pathway comparison
# ------------------------------------------------------------

rank_data <- rankNet(cellchat_merged, mode = "comparison", stacked = TRUE,
                     do.stat = TRUE, return.data = TRUE)

write.csv(rank_data$signaling.contribution,
          "results/TME/cellchat_subtype_pathway_rank_table.csv",
          row.names = FALSE)

pdf("results/TME/cellchat_subtype_pathway_rank_FULL.pdf", width = 9, height = 14)
print(rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE))
dev.off()


# ------------------------------------------------------------
# 7. PROPERLY interpret rank_data — sort by actual contribution,
#    not the raw unsorted table (head() alone is misleading,
#    since the raw table isn't ordered by importance)
# ------------------------------------------------------------

rc <- rank_data$signaling.contribution

cat("\n=== TOP 15 pathways by contribution (descending) ===\n")
print(rc %>% arrange(desc(contribution)) %>% head(15))

cat("\n=== Significant (p<0.05) pathways with contribution > 0, sorted ===\n")
sig_pathways <- rc %>%
  filter(pvalues < 0.05, contribution > 0) %>%
  arrange(desc(contribution))
print(sig_pathways)

write.csv(sig_pathways, "results/TME/cellchat_significant_pathways_sorted.csv",
          row.names = FALSE)

cat("\n=== Candidate pathways (TGFb/SPP1/CXCL/CCL/MIF) across subtypes ===\n")
candidate_table <- rc %>% filter(name %in% candidate_pathways) %>%
  arrange(name, group)
print(candidate_table)

write.csv(candidate_table, "results/TME/cellchat_candidate_pathways_table.csv",
          row.names = FALSE)


# ------------------------------------------------------------
# 8. Bubble plot — TME senders -> malignant states, per subtype
# ------------------------------------------------------------

common_groups <- Reduce(intersect, lapply(object.list, function(x) levels(x@idents)))
malignant_states <- intersect(
  c("Basal-like", "Luminal Progenitor-like", "Mature Luminal-like",
    "EMT/Mesenchymal", "Cycling/Proliferative"),
  common_groups
)
non_malignant <- setdiff(common_groups, c(malignant_states, "Normal_epithelial"))

pdf("results/TME/cellchat_subtype_bubble_comparison.pdf", width = 12, height = 10)
print(netVisual_bubble(
  cellchat_merged,
  sources.use = non_malignant,
  targets.use = malignant_states,
  comparison = c(1, 2, 3),
  angle.x = 45
))
dev.off()


cat("\n============================================\n")
cat("Key result files:\n")
cat("  - cellchat_significant_pathways_sorted.csv   <- properly sorted,\n")
cat("      significant pathways only. TOP of this file = your\n")
cat("      strongest candidate for the 'mechanistic axis'.\n")
cat("  - cellchat_candidate_pathways_table.csv      <- TGFb/SPP1/CXCL/\n")
cat("      CCL/MIF values across TNBC/HER2+/ER+, side by side.\n")
cat("  - cellchat_subtype_bubble_comparison.pdf\n")
cat("  - cellchat_subtype_pathway_rank_FULL.pdf\n")
cat("============================================\n")


# ============================================================
# [REVISION] 9. Midkine (MK) checks that do not require re-running CellChat
#     (peer-review addition)
#   A. MK communication probabilities stored in the three subtype objects:
#      full network and network restricted to plausibly annotated groups
#      (the full network contains groups with implausible SingleR labels,
#      e.g. iPS cells, neurons, chondrocytes and MSC)
#   B. patient-level, expression-based MK proxy
#      (MDK in stromal senders x mean of nine MK receptor genes in malignant cells)
# Output: results/TME/MK_flow_full_network.csv, MK_flow_restricted.csv,
#         MK_per_patient_restricted.csv, MK_per_patient_preliminary.csv,
#         MK_leave_one_patient_out_restricted.csv, MK_per_patient.pdf
# ============================================================

## ---------------- A. MK flow from the stored CellChat objects ----------------
stromal_senders <- c("Fibroblasts", "Endothelial_cells", "Smooth_muscle_cells")
mal_states      <- c("Basal-like", "Luminal Progenitor-like", "Mature Luminal-like",
                     "EMT/Mesenchymal", "Cycling/Proliferative")
keep_senders <- c(stromal_senders, mal_states, "Normal_epithelial")
keep_receivers <- mal_states

full_out <- list(); restricted_out <- list()
for (s in c("TNBC", "HER2", "ER")) {
  cc <- readRDS(sprintf("results/TME/cellchat_%s.rds", s))
  pr <- cc@netP$prob
  mk_all <- pr[, , "MK"]
  cat("\n==", s, "(full network) ==\nMK information flow:", signif(sum(mk_all), 3),
      "| share of all pathways:", round(sum(mk_all) / sum(pr), 4),
      "| non-zero edges:", sum(mk_all > 0), "\n")
  top <- as.data.frame(as.table(mk_all)); names(top) <- c("sender", "receiver", "prob")
  print(head(top %>% arrange(desc(prob)), 6), row.names = FALSE)
  full_out[[s]] <- data.frame(subtype = s, MK_flow = sum(mk_all),
                              share_of_all_pathways = sum(mk_all) / sum(pr),
                              nonzero_edges = sum(mk_all > 0))

  sn <- intersect(keep_senders, dimnames(pr)[[1]])
  rc <- intersect(keep_receivers, dimnames(pr)[[2]])
  mk_r <- pr[sn, rc, "MK"]; all_r <- pr[sn, rc, , drop = FALSE]
  restricted_out[[s]] <- data.frame(subtype = s, MK_flow = sum(mk_r),
                                    share_of_all_pathways = sum(mk_r) / sum(all_r),
                                    nonzero_edges = sum(mk_r > 0),
                                    n_senders = length(sn), n_receivers = length(rc))
  rm(cc); gc()
}
full_out <- do.call(rbind, full_out)
full_out$relative_flow <- round(full_out$MK_flow / sum(full_out$MK_flow), 3)
restricted_out <- do.call(rbind, restricted_out)
restricted_out$relative_flow <- round(restricted_out$MK_flow / sum(restricted_out$MK_flow), 3)
print(full_out, row.names = FALSE); print(restricted_out, row.names = FALSE)
write.csv(full_out,       "results/TME/MK_flow_full_network.csv", row.names = FALSE)
write.csv(restricted_out, "results/TME/MK_flow_restricted.csv",   row.names = FALSE)

## ---------------- B. patient-level MK proxy ----------------
obj_mk <- readRDS("results/TME/full_annotated_for_cellchat.rds")
wu_meta <- read.csv("data/GSE176078/GSE176078/metadata.csv")
sample_subtype_mk <- unique(wu_meta[, c("orig.ident", "subtype")])
obj_mk$subtype <- sample_subtype_mk$subtype[match(obj_mk$orig.ident, sample_subtype_mk$orig.ident)]
lab <- as.character(obj_mk$final_celltype)

mk_receptors <- c("SDC1", "SDC2", "SDC4", "PTPRZ1", "LRP1", "NCL", "ITGA4", "ITGA6", "ITGB1")
expr_mk <- GetAssayData(obj_mk, assay = "RNA", layer = "data")
mk_receptors <- intersect(mk_receptors, rownames(expr_mk))
stopifnot("MDK" %in% rownames(expr_mk))
mdk <- as.numeric(expr_mk["MDK", ])

patient_level <- function(senders, min_mal, min_snd) {
  res <- do.call(rbind, lapply(unique(obj_mk$orig.ident), function(p) {
    ip <- obj_mk$orig.ident == p
    im <- which(ip & lab %in% mal_states); is <- which(ip & lab %in% senders)
    data.frame(patient = p, subtype = unique(obj_mk$subtype[ip])[1],
               n_malignant = length(im), n_senders = length(is),
               MDK_senders = if (length(is) >= min_snd) mean(mdk[is]) else NA,
               receptor_malignant = if (length(im) >= min_mal)
                 mean(Matrix::rowMeans(expr_mk[mk_receptors, im, drop = FALSE])) else NA)
  }))
  res$MK_potential <- res$MDK_senders * res$receptor_malignant
  res$subtype <- factor(res$subtype, levels = c("TNBC", "HER2+", "ER+"))
  res
}

report_patient_level <- function(res) {
  ok <- res[!is.na(res$MK_potential), ]
  print(res %>% arrange(subtype, desc(MK_potential)), row.names = FALSE)
  cat("\nPatients analysed per subtype:\n"); print(table(ok$subtype))
  for (v in c("MK_potential", "MDK_senders", "receptor_malignant")) {
    cat("\n---", v, "| Kruskal-Wallis p =", signif(kruskal.test(ok[[v]] ~ ok$subtype)$p.value, 3), "\n")
    print(ok %>% group_by(subtype) %>%
            summarise(n = n(), median = round(median(.data[[v]]), 4),
                      min = round(min(.data[[v]]), 4), max = round(max(.data[[v]]), 4)))
  }
  cat("TNBC vs ER+ Wilcoxon p =", signif(wilcox.test(MK_potential ~ subtype,
      data = droplevels(ok[ok$subtype %in% c("TNBC", "ER+"), ]), exact = FALSE)$p.value, 3), "\n")
  loo <- do.call(rbind, lapply(seq_len(nrow(ok)), function(i) {
    m <- tapply(ok$MK_potential[-i], ok$subtype[-i], median)
    data.frame(dropped = ok$patient[i], TNBC = round(m["TNBC"], 3), HER2 = round(m["HER2+"], 3),
               ER = round(m["ER+"], 3), highest = names(which.max(m)))
  }))
  print(table(loo$highest))
  loo
}

# (i) preliminary analysis: senders = all non-malignant groups with MDK detected in >= 20% of cells
det_rate <- tapply(mdk > 0, lab, mean)
senders_auto <- setdiff(names(det_rate)[det_rate >= 0.20], c(mal_states, "Normal_epithelial"))
res_prelim <- patient_level(senders_auto, min_mal = 20, min_snd = 20)
invisible(report_patient_level(res_prelim))
write.csv(res_prelim, "results/TME/MK_per_patient_preliminary.csv", row.names = FALSE)

# (ii) main analysis: stromal senders shown in Fig. 5d; >= 50 malignant and >= 20 sender cells
res_mk <- patient_level(stromal_senders, min_mal = 50, min_snd = 20)
loo_mk <- report_patient_level(res_mk)
write.csv(res_mk, "results/TME/MK_per_patient_restricted.csv", row.names = FALSE)
write.csv(loo_mk, "results/TME/MK_leave_one_patient_out_restricted.csv", row.names = FALSE)

library(ggplot2)
ok_mk <- res_mk[!is.na(res_mk$MK_potential), ]
pdf("results/TME/MK_per_patient.pdf", width = 6, height = 5)
print(ggplot(ok_mk, aes(subtype, MK_potential)) +
        geom_boxplot(outlier.shape = NA, fill = NA) +
        geom_point(size = 2, position = position_jitter(width = 0.1, seed = 1)) +
        geom_text(aes(label = patient), size = 2.3, nudge_x = 0.2, check_overlap = TRUE) +
        labs(x = NULL, y = "MDK (stromal senders) x receptors (malignant)",
             title = "Per-patient MK potential") +
        theme_classic())
dev.off()
