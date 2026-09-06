# ============================================================
# FIGURE 4B — MALIGNANT STATE MODULE-SCORE HEATMAP
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

scores <- read_csv(
  "results/TME/malignant_state_scores_per_cell.csv",
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 2. Exact score columns from YOUR file
# ------------------------------------------------------------

score_cols <- c(
  "Basal_score",
  "LumProg_score",
  "MatureLum_score",
  "EMT_score",
  "Cycling_score"
)

# Check that all columns exist
missing_cols <- setdiff(score_cols, colnames(scores))

if (length(missing_cols) > 0) {
  stop(
    paste(
      "These columns are missing:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# 3. Calculate mean module score for each malignant state
# ------------------------------------------------------------

mean_scores <- scores %>%
  summarise(
    Basal_score      = mean(Basal_score, na.rm = TRUE),
    LumProg_score    = mean(LumProg_score, na.rm = TRUE),
    MatureLum_score  = mean(MatureLum_score, na.rm = TRUE),
    EMT_score        = mean(EMT_score, na.rm = TRUE),
    Cycling_score    = mean(Cycling_score, na.rm = TRUE)
  )

print(mean_scores)

# ------------------------------------------------------------
# 4. Convert to long format
# ------------------------------------------------------------

heatmap_df <- mean_scores %>%
  pivot_longer(
    cols = everything(),
    names_to = "State",
    values_to = "MeanScore"
  )

# Human-readable labels
heatmap_df$State <- factor(
  heatmap_df$State,
  levels = score_cols,
  labels = c(
    "Basal-like",
    "Luminal progenitor-like",
    "Mature luminal-like",
    "EMT/mesenchymal",
    "Cycling/proliferative"
  )
)

# ------------------------------------------------------------
# 5. IMPORTANT:
#    We need a 5 x 5 heatmap if we want
#    malignant state x module score.
#
#    So calculate average score FOR EACH malignant state
# ------------------------------------------------------------

state_heatmap <- scores %>%
  group_by(malignant_state) %>%
  summarise(
    Basal_score     = mean(Basal_score, na.rm = TRUE),
    LumProg_score   = mean(LumProg_score, na.rm = TRUE),
    MatureLum_score = mean(MatureLum_score, na.rm = TRUE),
    EMT_score       = mean(EMT_score, na.rm = TRUE),
    Cycling_score   = mean(Cycling_score, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 6. Convert to long format
# ------------------------------------------------------------

state_heatmap_long <- state_heatmap %>%
  pivot_longer(
    cols = all_of(score_cols),
    names_to = "Module",
    values_to = "Score"
  )

# ------------------------------------------------------------
# 7. Clean labels
# ------------------------------------------------------------

state_heatmap_long$Module <- factor(
  state_heatmap_long$Module,
  levels = score_cols,
  labels = c(
    "Basal-like",
    "Luminal progenitor-like",
    "Mature luminal-like",
    "EMT/mesenchymal",
    "Cycling/proliferative"
  )
)

# Clean malignant-state labels
state_heatmap_long$malignant_state <- factor(
  state_heatmap_long$malignant_state,
  levels = c(
    "Basal-like",
    "Luminal Progenitor-like",
    "Mature Luminal-like",
    "EMT/Mesenchymal",
    "Cycling/Proliferative"
  ),
  labels = c(
    "Basal-like",
    "Luminal progenitor-like",
    "Mature luminal-like",
    "EMT/mesenchymal",
    "Cycling/proliferative"
  )
)

# ------------------------------------------------------------
# 8. Plot heatmap
# ------------------------------------------------------------

p4B <- ggplot(
  state_heatmap_long,
  aes(
    x = Module,
    y = malignant_state,
    fill = Score
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.6
  ) +
  
  geom_text(
    aes(label = sprintf("%.2f", Score)),
    size = 3.5
  ) +
  
  scale_fill_gradient2(
    midpoint = 0,
    name = "Mean\nmodule score"
  ) +
  
  labs(
    title = "Malignant-state module scores",
    x = "Malignant transcriptional program",
    y = "Assigned malignant state"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

print(p4B)

# ------------------------------------------------------------
# 9. Save figure
# ------------------------------------------------------------

dir.create(
  "figures/Malignant",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  "figures/Malignant/Figure4B_Malignant_State_ModuleScore_Heatmap.png",
  p4B,
  width = 8,
  height = 6,
  dpi = 600
)

cat("\nFigure 4B generated successfully.\n")