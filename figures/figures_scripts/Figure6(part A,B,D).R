# ============================================================
# FIGURE 6A —60-GENE SIGNATURE SCHEMATIC
# ============================================================

library(ggplot2)
library(grid)

# ------------------------------------------------------------
# Create publication-style schematic
# ------------------------------------------------------------

p6A <- ggplot() +
  
  # ==========================================================
# LEFT: SINGLE-CELL DISCOVERY
# ==========================================================

annotate(
  "text",
  x = 1.6,
  y = 4.65,
  label = "Single-cell discovery",
  fontface = "bold",
  size = 6
) +
  
  # Cycling state
  annotate(
    "label",
    x = 1.6,
    y = 3.55,
    label =
      "Cycling / proliferative state\n\n50 marker genes\nMKI67 • TOP2A • CDK1 •\nPLK1 • UBE2C • AURKB",
    size = 4.2,
    lineheight = 1.15,
    label.padding = unit(0.8, "lines"),
    label.r = unit(0.15, "lines")
  ) +
  
  # MK pathway
  annotate(
    "label",
    x = 1.6,
    y = 1.55,
    label =
      "Midkine (MK) signaling\n\n10 pathway genes\nMDK • LRP1 • PTPRZ1 •\nSDC1 • SDC2 • SDC4",
    size = 4.2,
    lineheight = 1.15,
    label.padding = unit(0.8, "lines"),
    label.r = unit(0.15, "lines")
  ) +
  
  # ==========================================================
# CENTER: INTEGRATION
# ==========================================================

annotate(
  "text",
  x = 5,
  y = 4.65,
  label = "Signature integration",
  fontface = "bold",
  size = 6
) +
  
  annotate(
    "label",
    x = 5,
    y = 2.55,
    label =
      "50 malignant-state markers\n+\n10 MK pathway genes\n\n────────────────\n\n60-gene integrated signature",
    size = 4.6,
    lineheight = 1.2,
    fontface = "bold",
    label.padding = unit(1.0, "lines"),
    label.r = unit(0.15, "lines")
  ) +
  
  # Arrows from left to center
  annotate(
    "segment",
    x = 2.55,
    xend = 3.65,
    y = 3.55,
    yend = 3.05,
    linewidth = 1.2,
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  
  annotate(
    "segment",
    x = 2.55,
    xend = 3.65,
    y = 1.55,
    yend = 2.05,
    linewidth = 1.2,
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  
  # ==========================================================
# RIGHT: CLINICAL VALIDATION
# ==========================================================

annotate(
  "text",
  x = 8.35,
  y = 4.65,
  label = "Independent validation",
  fontface = "bold",
  size = 6
) +
  
  annotate(
    "label",
    x = 8.35,
    y = 2.55,
    label =
      "METABRIC cohort\n\nn = 1,980 patients\n\nSignature score\n↓\nOverall survival",
    size = 4.8,
    lineheight = 1.2,
    label.padding = unit(1.0, "lines"),
    label.r = unit(0.15, "lines")
  ) +
  
  # Arrow center → right
  annotate(
    "segment",
    x = 6.35,
    xend = 7.05,
    y = 2.55,
    yend = 2.55,
    linewidth = 1.2,
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  
  # ==========================================================
# FINAL EMPHASIS
# ==========================================================

annotate(
  "text",
  x = 5,
  y = 0.35,
  label =
    "Single-cell-derived signature translated to clinical prognosis",
  fontface = "italic",
  size = 4.3
) +
  
  theme_void() +
  
  xlim(0.2, 9.8) +
  ylim(0, 5.2)

# ------------------------------------------------------------
# Display
# ------------------------------------------------------------

print(p6A)

# ------------------------------------------------------------
# Save high-resolution outputs
# ------------------------------------------------------------

dir.create(
  "figures/Signature",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  "figures/Signature/Figure6A_60Gene_Signature_Construction.png",
  p6A,
  width = 13,
  height = 6,
  dpi = 600,
  bg = "white"
)


cat(
  "\nFigure 6A publication-style schematic generated successfully.\n"
)
# ============================================================
# FIGURE 6B
# METABRIC SIGNATURE-SCORE DISTRIBUTION
# ============================================================

bulk <- read_csv(
  "results/TME/bulk_signature_scores.csv",
  show_col_types = FALSE
)

cat("\n========================================\n")
cat("METABRIC signature-score file\n")
cat("========================================\n")

print(colnames(bulk))
print(head(bulk))

# ------------------------------------------------------------
# Identify numeric score column
# ------------------------------------------------------------

numeric_cols <- names(
  bulk[
    sapply(bulk, is.numeric)
  ]
)

# Remove obvious clinical variables if present
numeric_cols <- numeric_cols[
  !grepl(
    "age|year|time|months|status|event|death",
    numeric_cols,
    ignore.case = TRUE
  )
]

if (length(numeric_cols) == 0) {
  stop(
    "Could not identify the signature-score column automatically."
  )
}

score_col <- numeric_cols[1]

cat("\nUsing signature-score column:", score_col, "\n")

# ------------------------------------------------------------
# Clean data
# ------------------------------------------------------------

bulk_plot <- bulk %>%
  filter(!is.na(.data[[score_col]]))

median_score <- median(
  bulk_plot[[score_col]],
  na.rm = TRUE
)

cat("Median signature score:", median_score, "\n")

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p6B <- ggplot(
  bulk_plot,
  aes(x = .data[[score_col]])
) +
  
  geom_histogram(
    bins = 50,
    color = "white"
  ) +
  
  geom_vline(
    xintercept = median_score,
    linetype = "dashed",
    linewidth = 0.9
  ) +
  
  annotate(
    "text",
    x = median_score,
    y = Inf,
    label = "Median split",
    vjust = 1.5,
    hjust = -0.05,
    size = 4
  ) +
  
  labs(
    title = "METABRIC 60-gene signature score",
    x = "Signature score",
    y = "Number of patients"
  ) +
  
  theme_classic(base_size = 13)

print(p6B)

ggsave(
  "figures/Signature/Figure6B_METABRIC_Signature_Score_Distribution.png",
  p6B,
  width = 7,
  height = 5.5,
  dpi = 600
)
# ============================================================
# FIGURE 6D — COX FOREST PLOT
# ============================================================

library(survival)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)

# ------------------------------------------------------------
# 1. Load METABRIC data
# ------------------------------------------------------------

bulk <- read.csv(
  "results/TME/bulk_signature_scores.csv",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 2. Prepare variables
# ------------------------------------------------------------

bulk <- bulk %>%
  mutate(
    signature_score = as.numeric(signature_score),
    OS_MONTHS = as.numeric(OS_MONTHS),
    
    OS_event = case_when(
      grepl("1:DECEASED", OS_STATUS) ~ 1,
      grepl("0:LIVING", OS_STATUS) ~ 0,
      TRUE ~ NA_real_
    ),
    
    AGE_AT_DIAGNOSIS =
      as.numeric(AGE_AT_DIAGNOSIS),
    
    CLAUDIN_SUBTYPE =
      as.factor(CLAUDIN_SUBTYPE)
  )

# ------------------------------------------------------------
# 3. Complete cases
# ------------------------------------------------------------

cox_data <- bulk %>%
  filter(
    !is.na(signature_score),
    !is.na(OS_MONTHS),
    !is.na(OS_event)
  )

# ------------------------------------------------------------
# 4. Unadjusted Cox model
# ------------------------------------------------------------

cox_unadj <- coxph(
  Surv(OS_MONTHS, OS_event) ~ signature_score,
  data = cox_data
)

unadj_hr <- exp(coef(cox_unadj)["signature_score"])
unadj_ci <- exp(confint(cox_unadj)["signature_score", ])
unadj_p  <- summary(cox_unadj)$coefficients[
  "signature_score",
  "Pr(>|z|)"
]

# ------------------------------------------------------------
# 5. Adjusted Cox model
# ------------------------------------------------------------

cox_adj <- coxph(
  Surv(OS_MONTHS, OS_event) ~
    signature_score +
    AGE_AT_DIAGNOSIS +
    CLAUDIN_SUBTYPE,
  data = cox_data
)

adj_hr <- exp(coef(cox_adj)["signature_score"])
adj_ci <- exp(confint(cox_adj)["signature_score", ])
adj_p  <- summary(cox_adj)$coefficients[
  "signature_score",
  "Pr(>|z|)"
]

# ------------------------------------------------------------
# 6. Create results table
# ------------------------------------------------------------

forest_df <- data.frame(
  Model = c(
    "Unadjusted",
    "Adjusted for age + subtype"
  ),
  
  HR = c(
    unadj_hr,
    adj_hr
  ),
  
  CI_low = c(
    unadj_ci[1],
    adj_ci[1]
  ),
  
  CI_high = c(
    unadj_ci[2],
    adj_ci[2]
  ),
  
  P = c(
    unadj_p,
    adj_p
  )
)

print(forest_df)

# ------------------------------------------------------------
# 7. Formatting
# ------------------------------------------------------------

forest_df <- forest_df %>%
  mutate(
    Model = factor(
      Model,
      levels = rev(Model)
    ),
    
    CI_text = sprintf(
      "%.2f–%.2f",
      CI_low,
      CI_high
    ),
    
    HR_text = sprintf(
      "%.2f",
      HR
    ),
    
    P_text = ifelse(
      P < 0.001,
      formatC(P, format = "e", digits = 2),
      sprintf("%.3f", P)
    )
  )

# ------------------------------------------------------------
# 8. Forest plot
# ------------------------------------------------------------

p6D_forest <- ggplot(
  forest_df,
  aes(
    x = HR,
    y = Model
  )
) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  geom_errorbarh(
    aes(
      xmin = CI_low,
      xmax = CI_high
    ),
    height = 0.18,
    linewidth = 1
  ) +
  
  geom_point(
    size = 4
  ) +
  
  scale_x_log10(
    breaks = c(0.75, 1, 1.25, 1.5, 2),
    limits = c(0.7, 2.1)
  ) +
  
  labs(
    title = "Cox regression",
    x = "Hazard ratio (log scale)",
    y = NULL
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text.y = element_text(
      face = "bold"
    )
  )

# ------------------------------------------------------------
# 9. Results table
# ------------------------------------------------------------

table_df <- forest_df %>%
  mutate(
    Model = as.character(Model)
  ) %>%
  select(
    Model,
    HR,
    CI_text,
    P_text
  )

table_grob <- tableGrob(
  table_df,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 11,
    core = list(
      fg_params = list(hjust = 0.5)
    ),
    colhead = list(
      fg_params = list(fontface = "bold")
    )
  )
)

# ------------------------------------------------------------
# 10. Combine forest plot + table
# ------------------------------------------------------------

combined_plot <- arrangeGrob(
  p6D_forest,
  table_grob,
  ncol = 1,
  heights = c(3.5, 1.4)
)

# ------------------------------------------------------------
# 11. Save
# ------------------------------------------------------------

dir.create(
  "figures/Signature",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  "figures/Signature/Figure6D_Publication_Cox_Plot.png",
  combined_plot,
  width = 9,
  height = 6.5,
  dpi = 600,
  bg = "white"
)

cat("\nFigure 6D publication-style plot generated successfully.\n")