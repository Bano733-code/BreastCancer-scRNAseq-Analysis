# ============================================================
# MANUSCRIPT FIGURE 3 (Identification of CNV-high epithelial cells)
#   Fig. 3a  epithelial cells on the original UMAP -> script 17, section 7
#            (figures/Malignant/Epithelial_Cells_Original_UMAP.png)
#   Fig. 3b  epithelial cells per sample (26 samples) -> script 17, section 8
#            (figures/Malignant/Epithelial_Cells_By_Sample.png)
#   Fig. 3c  distribution of chromosome-level CNV scores -> section at the END
#            of this file (figures/Malignant/Figure3C_CNV_score_distribution.png)
# The two plots below (CNV-high proportion per patient and CNV-high vs CNV-low
# classification) are supporting figures that are not part of the final manuscript figure.
# Run from the project root (no setwd() is used).
# ============================================================

# ============================================================
# FIGURE 3B — PATIENT-LEVEL CNV-HIGH PROPORTION
# ============================================================

library(readr)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load CNV results
# ------------------------------------------------------------

cnv <- read_csv(
  "results/TME/cnv_scores_per_cell.csv",
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 2. Create CNV-high indicator
# ------------------------------------------------------------

cnv <- cnv %>%
  mutate(
    CNV_high = ifelse(
      cnv_call == "CNV-high (malignant candidate)",
      1,
      0
    )
  )

# ------------------------------------------------------------
# 3. Calculate CNV-high proportion for every sample
# ------------------------------------------------------------

sample_cnv <- cnv %>%
  group_by(sample) %>%
  summarise(
    total_epithelial = n(),
    CNV_high_cells = sum(CNV_high),
    CNV_high_percent =
      100 * CNV_high_cells / total_epithelial,
    .groups = "drop"
  ) %>%
  arrange(desc(CNV_high_percent))

# ------------------------------------------------------------
# 4. Print results
# ------------------------------------------------------------

cat("\nPatient-level CNV-high proportions:\n")
print(sample_cnv)

cat(
  "\nMinimum CNV-high proportion:",
  round(min(sample_cnv$CNV_high_percent), 2),
  "%\n"
)

cat(
  "Maximum CNV-high proportion:",
  round(max(sample_cnv$CNV_high_percent), 2),
  "%\n"
)

# ------------------------------------------------------------
# 5. Make sample order follow percentage
# ------------------------------------------------------------

sample_cnv$sample <- factor(
  sample_cnv$sample,
  levels = sample_cnv$sample
)

# ------------------------------------------------------------
# 6. Publication-style horizontal bar plot
# ------------------------------------------------------------

p3B <- ggplot(
  sample_cnv,
  aes(
    x = sample,
    y = CNV_high_percent
  )
) +

  geom_col(
    width = 0.7
  ) +

  coord_flip() +

  geom_text(
    aes(
      label = paste0(
        round(CNV_high_percent, 1),
        "%"
      )
    ),
    hjust = -0.15,
    size = 3
  ) +

  scale_y_continuous(
    limits = c(
      0,
      max(sample_cnv$CNV_high_percent) * 1.12
    ),
    expand = c(0, 0)
  ) +

  labs(
    title = "CNV-high epithelial cells across patient samples",
    x = "Patient sample",
    y = "CNV-high epithelial cells (%)"
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 13
    ),
    axis.text.y = element_text(
      size = 8
    ),
    axis.title = element_text(
      face = "bold"
    )
  )

print(p3B)

# ------------------------------------------------------------
# 7. Save high-resolution PNG
# ------------------------------------------------------------

dir.create(
  "figures/Malignant",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  "figures/Malignant/Figure3B_CNV_High_Proportion_By_Sample.png",
  p3B,
  width = 7,
  height = 8,
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# 8. Save the summary table
# ------------------------------------------------------------

write.csv(
  sample_cnv,
  "results/TME/CNV_High_Proportion_By_Sample.csv",
  row.names = FALSE
)

cat(
  "\nFigure 3B generated successfully.\n"
)
# ============================================================
# FIGURE 3C — CNV-HIGH VS CNV-LOW CLASSIFICATION
# ============================================================

# ------------------------------------------------------------
# Count cells in each CNV category
# ------------------------------------------------------------

cnv_counts <- cnv %>%
  count(CNV_group) %>%
  mutate(
    Percentage = 100 * n / sum(n)
  )

print(cnv_counts)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p3C <- ggplot(
  cnv_counts,
  aes(x = CNV_group, y = n)
) +
  
  geom_col(
    width = 0.65,
    fill = "grey65"
  ) +
  
  geom_text(
    aes(
      label = paste0(
        format(n, big.mark = ","),
        "\n(",
        round(Percentage, 1),
        "%)"
      )
    ),
    vjust = -0.25,
    size = 4.5
  ) +
  
  labs(
    title = "Classification of epithelial cells by CNV score",
    x = NULL,
    y = "Number of epithelial cells"
  ) +
  
  expand_limits(y = max(cnv_counts$n) * 1.15) +
  
  theme_classic(base_size = 13)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

ggsave(
  "figures/Malignant/Figure3C_CNV_High_vs_Low.png",
  p3C,
  width = 6.5,
  height = 5.5,
  dpi = 600
)


# ============================================================
# [REVISION] MANUSCRIPT FIGURE 3c (final): distribution of chromosome-level CNV scores
# The title now says "Chromosome-Level" (the score is the SD across chromosomes of
# the mean reference-adjusted expression; chromosome arms are not resolved).
# The y-axis is the probability density (area under the curve = 1), not a cell count.
# ============================================================

library(readr)
library(dplyr)
library(ggplot2)

dir.create("figures/Malignant", recursive = TRUE, showWarnings = FALSE)

cnv_scores <- read_csv("results/TME/cnv_scores_per_cell.csv", show_col_types = FALSE)
threshold  <- 0.025   # empirical 95th percentile of the reference-cell scores (script 18)

p3c <- ggplot(cnv_scores, aes(x = cnv_score)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80,
                 fill = "grey75", colour = "white") +
  geom_density(linewidth = 1) +
  geom_vline(xintercept = threshold, linetype = "dashed", linewidth = 1) +
  annotate("text", x = threshold + 0.001, y = Inf, label = "Threshold = 0.025",
           hjust = 0, vjust = 1.6, size = 5) +
  labs(title = "Distribution of Chromosome-Level CNV Scores",
       x = "CNV Score", y = "Density") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))

ggsave("figures/Malignant/Figure3C_CNV_score_distribution.png", p3c,
       width = 8, height = 6, dpi = 300, bg = "white")
print(p3c)
