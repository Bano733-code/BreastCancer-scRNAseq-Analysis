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
