# ============================================================
# Survival Analysis: Kaplan-Meier + Cox regression
# ============================================================

# install.packages(c("survival", "survminer"))  # run once if needed
library(survival)
library(survminer)
library(dplyr)


# ------------------------------------------------------------
# 1. Load merged score + clinical data from Step 9
# ------------------------------------------------------------

merged <- read.csv("results/TME/bulk_signature_scores.csv")

cat("Columns available:\n")
print(colnames(merged))


# ------------------------------------------------------------
# 2. Identify and clean survival columns
#
#    METABRIC clinical files typically have:
#      OS_MONTHS   (overall survival time in months)
#      OS_STATUS   (e.g. "0:LIVING" / "1:DECEASED")
#
#    TCGA files may use similar or slightly different naming
#    (check your actual column names from the print() above
#     and adjust the lines below if they differ)
# ------------------------------------------------------------

stopifnot("OS_MONTHS" %in% colnames(merged))
stopifnot("OS_STATUS" %in% colnames(merged))

merged <- merged %>%
  filter(!is.na(OS_MONTHS), !is.na(OS_STATUS))

# Convert OS_STATUS to a numeric event indicator (1 = death, 0 = alive)
merged$event <- ifelse(grepl("DECEASED|1:", merged$OS_STATUS), 1, 0)

cat("\nEvent distribution:\n")
print(table(merged$event))

cat("\nFollow-up time summary (months):\n")
print(summary(merged$OS_MONTHS))


# ------------------------------------------------------------
# 3. Kaplan-Meier: High vs Low signature score
# ------------------------------------------------------------

merged$group <- factor(merged$group, levels = c("Low", "High"))

surv_obj <- Surv(time = merged$OS_MONTHS, event = merged$event)
km_fit <- survfit(surv_obj ~ group, data = merged)

# Log-rank test p-value
logrank_p <- survdiff(surv_obj ~ group, data = merged)
cat("\nLog-rank test:\n")
print(logrank_p)

# Plot
km_plot <- ggsurvplot(
  km_fit,
  data = merged,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  xlab = "Time (months)",
  ylab = "Overall Survival Probability",
  title = "Gene Signature: High vs Low — Overall Survival",
  legend.title = "Signature",
  legend.labs = c("Low", "High"),
  palette = c("#2C7BB6", "#D7191C")
)

pdf("results/TME/survival_KM_signature.pdf", width = 8, height = 8)
print(km_plot)
dev.off()

cat("\nSaved: results/TME/survival_KM_signature.pdf\n")


# ------------------------------------------------------------
# 4. Cox regression — UNADJUSTED (signature alone)
# ------------------------------------------------------------

cox_unadj <- coxph(surv_obj ~ group, data = merged)
cat("\n=== Unadjusted Cox regression ===\n")
print(summary(cox_unadj))


# ------------------------------------------------------------
# 5. Cox regression — ADJUSTED for stage/subtype (if available)
#
#    Check merged column names above for the exact stage/subtype
#    column your clinical file uses (commonly one of:
#    TUMOR_STAGE, AJCC_PATHOLOGIC_TUMOR_STAGE, CLAUDIN_SUBTYPE,
#    or the "subtype" column you may have merged in yourself)
# ------------------------------------------------------------

adjustment_vars <- intersect(
  c("AGE_AT_DIAGNOSIS", "TUMOR_STAGE", "CLAUDIN_SUBTYPE",
    "AJCC_PATHOLOGIC_TUMOR_STAGE"),
  colnames(merged)
)

cat("\nAdjustment variables found and used:", paste(adjustment_vars, collapse = ", "), "\n")

if (length(adjustment_vars) > 0) {
  
  formula_str <- paste("surv_obj ~ group +", paste(adjustment_vars, collapse = " + "))
  cox_adj <- coxph(as.formula(formula_str), data = merged)
  
  cat("\n=== Adjusted Cox regression ===\n")
  print(summary(cox_adj))
  
  cox_summary <- data.frame(
    model = c("Unadjusted", "Adjusted"),
    HR = c(exp(coef(cox_unadj))["groupHigh"], exp(coef(cox_adj))["groupHigh"]),
    p_value = c(summary(cox_unadj)$coefficients["groupHigh", "Pr(>|z|)"],
                summary(cox_adj)$coefficients["groupHigh", "Pr(>|z|)"])
  )
  
} else {
  cat("\nNo standard stage/subtype columns found — check your clinical\n")
  cat("file's column names and update 'adjustment_vars' above manually.\n")
  
  cox_summary <- data.frame(
    model = "Unadjusted",
    HR = exp(coef(cox_unadj))["groupHigh"],
    p_value = summary(cox_unadj)$coefficients["groupHigh", "Pr(>|z|)"]
  )
}

print(cox_summary)
write.csv(cox_summary, "results/TME/survival_cox_summary.csv", row.names = FALSE)


# ------------------------------------------------------------
# 6. Continuous score version (robustness check — treats
#    signature score as continuous rather than High/Low binary,
#    avoids any concern about median-split cutpoint sensitivity)
# ------------------------------------------------------------

cox_continuous <- coxph(surv_obj ~ signature_score, data = merged)
cat("\n=== Cox regression, CONTINUOUS signature score ===\n")
print(summary(cox_continuous))


cat("\n============================================\n")
cat("Saved: results/TME/survival_KM_signature.pdf\n")
cat("Saved: results/TME/survival_cox_summary.csv\n")
cat("============================================\n")
cat("\nPROJECT COMPLETE — all steps done.\n")
