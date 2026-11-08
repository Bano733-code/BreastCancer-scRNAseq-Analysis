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

# [REVISION] NOTE: count() returns only (sample, state) combinations that
# exist. A sample that has no cell of a given state therefore does NOT enter
# the test for that state in the primary analysis below (2-8 samples per
# state). The sensitivity analysis at the end of this script keeps all
# patients with CNV-high cells and sets absent states to zero.
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

if (!requireNamespace("FSA", quietly = TRUE)) install.packages("FSA")   # [REVISION] install only if missing
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


# ============================================================
# [REVISION] 10. Sensitivity analysis: all patients kept (absent state = 0)
#     (peer-review addition). Adds Kruskal-Wallis (BH), Dunn post-hoc tests,
#     leave-one-patient-out analysis, per-patient medians and a per-patient
#     Shannon diversity index (quantifies "broadest distribution of states").
# Output: results/TME/malignant_state_per_patient.csv, S7_*.csv,
#         results/TME/malignant_state_shannon.csv
# ============================================================

library(tidyr)
if (!requireNamespace("FSA", quietly = TRUE)) install.packages("FSA")

md_states <- malignant_obj@meta.data[, c("orig.ident", "subtype", "malignant_state")]
pp <- md_states %>%
  count(orig.ident, subtype, malignant_state) %>%
  group_by(orig.ident, subtype) %>%
  mutate(n_pat = sum(n), prop = n / n_pat) %>%
  ungroup()
wide <- pp %>%
  select(orig.ident, subtype, n_pat, malignant_state, prop) %>%
  pivot_wider(names_from = malignant_state, values_from = prop, values_fill = 0)
states <- setdiff(names(wide), c("orig.ident", "subtype", "n_pat"))
cat("\nPatients per subtype with CNV-high cells:\n"); print(table(wide$subtype))   # ER+ 8, HER2+ 4, TNBC 7

# 10a. Kruskal-Wallis per state (BH across the five states)
kw <- sapply(states, function(s) kruskal.test(wide[[s]] ~ factor(wide$subtype))$p.value)
kw_tab <- data.frame(state = states, KW_p = signif(kw, 3), BH_adj_p = signif(p.adjust(kw, "BH"), 3))
print(kw_tab)

# 10b. Dunn post-hoc tests (BH)
dunn_all <- do.call(rbind, lapply(states, function(s)
  cbind(state = s, FSA::dunnTest(wide[[s]] ~ factor(wide$subtype), method = "bh")$res)))
print(dunn_all, row.names = FALSE)

# 10c. Leave-one-patient-out (BH across states)
loo <- do.call(rbind, lapply(seq_len(nrow(wide)), function(i) {
  w <- wide[-i, ]
  p <- sapply(states, function(s) kruskal.test(w[[s]] ~ factor(w$subtype))$p.value)
  a <- p.adjust(p, "BH")
  data.frame(dropped = wide$orig.ident[i], n_sig = sum(a < 0.05),
             max_adj_p = round(max(a), 3), min_adj_p = round(min(a), 3))
}))
print(loo, row.names = FALSE)
cat("Leave-one-out: all five states significant in", sum(loo$n_sig == 5), "of", nrow(loo),
    "runs; range of n_sig:", range(loo$n_sig), "\n")

# 10d. Per-patient medians
print(as.data.frame(wide %>% group_by(subtype) %>%
        summarise(across(all_of(states), ~ round(median(.x), 3)), n_pat = n())))

# 10e. Shannon diversity per patient (patients with >= 20 malignant cells)
wd <- wide %>% filter(n_pat >= 20)
wd$shannon      <- apply(wd[, states], 1, function(p) { p <- p[p > 0]; -sum(p * log(p)) })
wd$shannon_norm <- wd$shannon / log(length(states))
print(as.data.frame(wd %>% group_by(subtype) %>%
        summarise(n = n(), median = round(median(shannon), 3),
                  min = round(min(shannon), 3), max = round(max(shannon), 3))))
cat("Kruskal-Wallis p (Shannon):", signif(kruskal.test(shannon ~ factor(subtype), data = wd)$p.value, 3), "\n")
print(FSA::dunnTest(shannon ~ factor(subtype), data = wd, method = "bh")$res[, c("Comparison", "Z", "P.adj")])

write.csv(wide,     "results/TME/malignant_state_per_patient.csv", row.names = FALSE)
write.csv(kw_tab,   "results/TME/S7_KW_zero_inclusive.csv",        row.names = FALSE)
write.csv(dunn_all, "results/TME/S7_Dunn_zero_inclusive.csv",      row.names = FALSE)
write.csv(loo,      "results/TME/S7_leave_one_out.csv",            row.names = FALSE)
write.csv(wd,       "results/TME/malignant_state_shannon.csv",     row.names = FALSE)
