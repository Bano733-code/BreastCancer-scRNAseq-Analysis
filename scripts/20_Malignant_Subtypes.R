# ============================================================
# Compare Malignant State Proportions Across Molecular Subtypes
# (TNBC / HER2+ / ER+)
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load malignant states object (Step 3 output)
# ------------------------------------------------------------

malignant_obj <- readRDS("results/TME/malignant_states_classified.rds")

cat("Malignant cells:", ncol(malignant_obj), "\n")
print(table(malignant_obj$malignant_state))


# ------------------------------------------------------------
# 2. Load original author-provided subtype metadata
#    and build a clean per-sample lookup table
# ------------------------------------------------------------

subtype_meta <- read.csv("data/GSE176078/GSE176078/metadata.csv")

sample_subtype <- unique(subtype_meta[, c("orig.ident", "subtype")])

cat("\nUnique samples with subtype:\n")
print(sample_subtype)

# Safety check: each sample should map to exactly ONE subtype
stopifnot(!any(duplicated(sample_subtype$orig.ident)))

cat("\nSample count per subtype:\n")
print(table(sample_subtype$subtype))


# ------------------------------------------------------------
# 3. Merge subtype into the malignant object
# ------------------------------------------------------------

malignant_obj$subtype <- sample_subtype$subtype[
  match(malignant_obj$orig.ident, sample_subtype$orig.ident)
]

# Confirm no missing subtype assignments
cat("\nCells with missing subtype (should be 0):",
    sum(is.na(malignant_obj$subtype)), "\n")

print(table(malignant_obj$orig.ident, malignant_obj$subtype))


# ------------------------------------------------------------
# 4. Malignant state proportions per subtype
# ------------------------------------------------------------

state_subtype_tab <- table(malignant_obj$subtype, malignant_obj$malignant_state)
print(state_subtype_tab)

# Convert to proportions (row-wise, i.e. within each subtype)
state_subtype_prop <- prop.table(state_subtype_tab, margin = 1)
print(round(state_subtype_prop * 100, 1))


# ------------------------------------------------------------
# 5. Stacked bar plot — malignant state composition per subtype
# ------------------------------------------------------------

plot_df <- as.data.frame(state_subtype_tab)
colnames(plot_df) <- c("subtype", "malignant_state", "n_cells")

ggplot(plot_df, aes(x = subtype, y = n_cells, fill = malignant_state)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Malignant Cell State Composition by Molecular Subtype",
    x = "Molecular Subtype",
    y = "Proportion of Malignant Cells",
    fill = "Malignant State"
  ) +
  theme_minimal(base_size = 13)

ggsave("results/TME/malignant_state_by_subtype_stackedbar.png",
       width = 7, height = 5, dpi = 300)


# ------------------------------------------------------------
# 6. Statistical test — is state distribution subtype-dependent?
#    Chi-square on the full contingency table first.
#    Falls back to Fisher's exact test if any expected cell
#    count is too low (<5), since chi-square becomes unreliable.
# ------------------------------------------------------------

chi_test <- chisq.test(state_subtype_tab)
print(chi_test)

low_expected <- any(chi_test$expected < 5)
cat("\nAny expected cell count < 5?:", low_expected, "\n")

if (low_expected) {
  cat("Running Fisher's exact test (simulated p-value, since table is large)...\n")
  fisher_test <- fisher.test(state_subtype_tab, simulate.p.value = TRUE, B = 10000)
  print(fisher_test)
}


# ------------------------------------------------------------
# 7. Pairwise subtype comparisons per malignant state
#    (which specific state drives the overall difference?)
# ------------------------------------------------------------

states <- unique(malignant_obj$malignant_state)
subtypes <- unique(malignant_obj$subtype)

pairwise_results <- list()

for (st in states) {
  # 2x3 table: this state vs "all other states", across subtypes
  is_state <- malignant_obj$malignant_state == st
  tab <- table(malignant_obj$subtype, is_state)
  
  test <- tryCatch(
    chisq.test(tab),
    warning = function(w) fisher.test(tab, simulate.p.value = TRUE, B = 10000)
  )
  
  pairwise_results[[st]] <- data.frame(
    state = st,
    p_value = test$p.value
  )
}

pairwise_df <- do.call(rbind, pairwise_results)
pairwise_df$p_adj <- p.adjust(pairwise_df$p_value, method = "BH")

cat("\nPer-state subtype association (BH-adjusted):\n")
print(pairwise_df)

write.csv(pairwise_df, "results/TME/malignant_state_subtype_association.csv",
          row.names = FALSE)


# ------------------------------------------------------------
# 8. Sample-level proportions (for a boxplot / more rigorous
#    stats that account for sample-level variability rather
#    than pooling all cells as independent — recommended for
#    the write-up, since cells within a sample aren't independent)
# ------------------------------------------------------------

sample_state_prop <- malignant_obj@meta.data %>%
  count(orig.ident, subtype, malignant_state) %>%
  group_by(orig.ident) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

write.csv(sample_state_prop, "results/TME/sample_level_state_proportions.csv",
          row.names = FALSE)

ggplot(sample_state_prop, aes(x = subtype, y = prop, fill = subtype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  facet_wrap(~ malignant_state, scales = "free_y") +
  labs(
    title = "Per-Sample Malignant State Proportions by Subtype",
    x = "Molecular Subtype", y = "Proportion of malignant cells (per sample)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave("results/TME/malignant_state_by_subtype_persample_boxplot.png",
       width = 9, height = 6, dpi = 300)


# ------------------------------------------------------------
# 9. Save final object + summary
# ------------------------------------------------------------

saveRDS(malignant_obj, "results/TME/malignant_states_with_subtype.rds")

cat("\n============================================\n")
cat("PHASE 1 STEP 4 COMPLETE\n")
cat("Saved: results/TME/malignant_states_with_subtype.rds\n")
cat("Saved: results/TME/malignant_state_by_subtype_stackedbar.png\n")
cat("Saved: results/TME/malignant_state_by_subtype_persample_boxplot.png\n")
cat("Saved: results/TME/malignant_state_subtype_association.csv\n")
cat("Saved: results/TME/sample_level_state_proportions.csv\n")
cat("============================================\n")

# ============================================================
# PHASE 1 — STEP 4 (ADDENDUM)
# Multiple-testing correction for the sample-level Kruskal-Wallis
# tests + Dunn's post-hoc pairwise comparisons
#
# Run this AFTER the main Step 4 script (needs sample_state_prop
# to already exist in your session)
# ============================================================

install.packages("FSA")   # run once if not already installed
library(FSA)
library(dplyr)


# ------------------------------------------------------------
# 1. Re-run Kruskal-Wallis per state, storing p-values properly
#    (instead of hand-typing them)
# ------------------------------------------------------------

kw_results <- list()

for (st in unique(sample_state_prop$malignant_state)) {
  sub_df <- sample_state_prop %>% filter(malignant_state == st)
  test <- kruskal.test(prop ~ subtype, data = sub_df)
  kw_results[[st]] <- data.frame(
    state   = st,
    chisq   = unname(test$statistic),
    df      = unname(test$parameter),
    p_value = test$p.value
  )
}

kw_df <- do.call(rbind, kw_results)
kw_df$p_adj <- p.adjust(kw_df$p_value, method = "BH")

cat("Sample-level Kruskal-Wallis results (BH-adjusted):\n")
print(kw_df)

write.csv(kw_df, "results/TME/malignant_state_kruskal_BHadjusted.csv",
          row.names = FALSE)


# ------------------------------------------------------------
# 2. Dunn's post-hoc test — WHICH subtype pairs differ,
#    for each malignant state (only meaningful for states
#    where the Kruskal-Wallis test above was significant)
# ------------------------------------------------------------

dunn_results <- list()

for (st in unique(sample_state_prop$malignant_state)) {
  
  sub_df <- sample_state_prop %>% filter(malignant_state == st)
  
  dunn_out <- dunnTest(prop ~ subtype, data = sub_df, method = "bh")
  
  dunn_table <- dunn_out$res
  dunn_table$state <- st
  
  dunn_results[[st]] <- dunn_table
}

dunn_df <- do.call(rbind, dunn_results)

cat("\nDunn's post-hoc pairwise comparisons (BH-adjusted within each state):\n")
print(dunn_df)

write.csv(dunn_df, "results/TME/malignant_state_dunn_posthoc.csv",
          row.names = FALSE)


# ------------------------------------------------------------
# 3. Quick summary — flag which pairwise comparisons are
#    significant after correction
# ------------------------------------------------------------

significant_pairs <- dunn_df %>% filter(P.adj < 0.05)

cat("\nSignificant pairwise subtype differences (adj. p < 0.05):\n")
print(significant_pairs[, c("state", "Comparison", "P.adj")])

cat("\n============================================\n")
cat("Saved: results/TME/malignant_state_kruskal_BHadjusted.csv\n")
cat("Saved: results/TME/malignant_state_dunn_posthoc.csv\n")
cat("============================================\n")
