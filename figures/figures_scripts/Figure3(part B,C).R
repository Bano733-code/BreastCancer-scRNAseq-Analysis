# ============================================================
# FIGURE 3B — CNV SCORE DISTRIBUTION
# ============================================================

library(ggplot2)
library(dplyr)

# ------------------------------------------------------------
# 1. Load CNV results
# ------------------------------------------------------------

cnv <- read.csv("results/TME/cnv_scores_per_cell.csv",
                stringsAsFactors = FALSE)

# Check
head(cnv)
table(cnv$cnv_call)

# ------------------------------------------------------------
# 2. Make clean CNV category
# ------------------------------------------------------------

cnv$CNV_group <- ifelse(
  cnv$cnv_call == "CNV-high (malignant candidate)",
  "CNV-high",
  "CNV-low"
)

# ------------------------------------------------------------
# 3. Plot distribution
# ------------------------------------------------------------

p3B <- ggplot(cnv, aes(x = cnv_score)) +
  
  geom_histogram(
    bins = 80,
    fill = "grey75",
    color = "white"
  ) +
  
  geom_density(
    aes(y = after_stat(count)),
    linewidth = 0.8
  ) +
  
  geom_vline(
    xintercept = 0.025,
    linetype = "dashed",
    linewidth = 0.9
  ) +
  
  annotate(
    "text",
    x = 0.025,
    y = Inf,
    label = "95th percentile threshold = 0.025",
    vjust = 1.5,
    hjust = -0.05,
    size = 4
  ) +
  
  labs(
    title = "Distribution of chromosome-arm CNV scores",
    x = "CNV score",
    y = "Number of epithelial cells"
  ) +
  
  theme_classic(base_size = 13)

# ------------------------------------------------------------
# 4. Save
# ------------------------------------------------------------

dir.create(
  "figures/Malignant",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  "figures/Malignant/Figure3B_CNV_Score_Distribution.png",
  p3B,
  width = 7,
  height = 5.5,
  dpi = 600
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