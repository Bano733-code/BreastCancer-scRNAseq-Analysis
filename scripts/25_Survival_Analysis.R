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

# [REVISION] Hazard ratios here are per 1-UNIT increase of the score; section 7
# (end of this script) also reports them per 1 SD (SD = 0.6626 units).
cox_continuous <- coxph(surv_obj ~ signature_score, data = merged)
cat("\n=== Cox regression, CONTINUOUS signature score ===\n")
print(summary(cox_continuous))


cat("\n============================================\n")
cat("Saved: results/TME/survival_KM_signature.pdf\n")
cat("Saved: results/TME/survival_cox_summary.csv\n")
cat("============================================\n")
cat("\nPROJECT COMPLETE — all steps done.\n")


# ============================================================
# [REVISION] 7. Extended survival analysis (peer-review addition)
#   - overall survival (OS), OS truncated at 120 months, disease-specific survival (DSS)
#   - median split, continuous score per 1 SD and per 1 unit
#   - three adjustment levels: unadjusted; age + subtype; full
#     (age, subtype, positive lymph nodes, chemotherapy, hormone therapy, radiotherapy)
#   - Schoenfeld-residual tests (cox.zph) for every model
# Output: results/TME/survival_models_full.csv, schoenfeld_full_model.pdf
# ============================================================

merged <- read.csv("results/TME/bulk_signature_scores.csv") %>%
  filter(!is.na(OS_MONTHS), !is.na(OS_STATUS))
merged$event_os <- ifelse(grepl("DECEASED|1:", merged$OS_STATUS), 1, 0)
merged$group    <- factor(merged$group, levels = c("Low", "High"))
merged$score_sd <- as.numeric(scale(merged$signature_score))             # per 1 SD
cat("SD of the signature score:", sd(merged$signature_score), "\n")      # 0.6626

# disease-specific survival: only "Died of Disease" counts as an event
merged$event_dss <- ifelse(merged$VITAL_STATUS == "Died of Disease", 1, 0)
merged$event_dss[is.na(merged$VITAL_STATUS)] <- NA
# OS truncated at 120 months (the Kaplan-Meier curves cross at about 250 months)
merged$time_120  <- pmin(merged$OS_MONTHS, 120)
merged$event_120 <- ifelse(merged$OS_MONTHS > 120, 0, merged$event_os)

for (v in intersect(c("CLAUDIN_SUBTYPE", "CHEMOTHERAPY", "HORMONE_THERAPY", "RADIO_THERAPY"), names(merged)))
  merged[[v]] <- factor(merged[[v]])

base_cov <- intersect(c("AGE_AT_DIAGNOSIS", "CLAUDIN_SUBTYPE"), names(merged))
full_cov <- intersect(c(base_cov, "LYMPH_NODES_EXAMINED_POSITIVE", "CHEMOTHERAPY",
                        "HORMONE_THERAPY", "RADIO_THERAPY"), names(merged))
cat("Base covariates:", base_cov, "\nFull covariates:", full_cov, "\n")
# Tumor grade and size are not available in this clinical export.

fit_one <- function(tm, ev, term, covs, label, endpoint) {
  d <- merged[complete.cases(merged[, c(tm, ev, term, covs)]), ]
  m <- coxph(as.formula(paste0("Surv(", tm, ",", ev, ") ~ ", paste(c(term, covs), collapse = " + "))), data = d)
  s <- summary(m); i <- grep(paste0("^", term), rownames(s$coefficients))[1]
  z <- cox.zph(m)
  data.frame(endpoint = endpoint, term = term, adjustment = label, N = nrow(d), events = sum(d[[ev]]),
             HR = round(s$coefficients[i, "exp(coef)"], 5),
             CI_low = round(s$conf.int[i, "lower .95"], 5), CI_high = round(s$conf.int[i, "upper .95"], 5),
             p = signif(s$coefficients[i, "Pr(>|z|)"], 3),
             PH_term_p = signif(z$table[term, "p"], 3), PH_global_p = signif(z$table["GLOBAL", "p"], 3))
}
endpoints <- list(list("OS_MONTHS", "event_os", "OS"),
                  list("time_120", "event_120", "OS_120m"),
                  list("OS_MONTHS", "event_dss", "DSS"))
cov_sets <- list(unadjusted = character(0), age_subtype = base_cov, full = full_cov)

models_full <- do.call(rbind, lapply(endpoints, function(e)
  do.call(rbind, lapply(c("group", "score_sd", "signature_score"), function(t)
    do.call(rbind, lapply(names(cov_sets), function(a)
      fit_one(e[[1]], e[[2]], t, cov_sets[[a]], a, e[[3]])))))))
print(models_full, row.names = FALSE, digits = 5)
write.csv(models_full, "results/TME/survival_models_full.csv", row.names = FALSE)

# Schoenfeld residuals, fully adjusted OS model (median split)
d_full <- merged[complete.cases(merged[, c("OS_MONTHS", "event_os", "group", full_cov)]), ]
m_full <- coxph(as.formula(paste("Surv(OS_MONTHS, event_os) ~ group +", paste(full_cov, collapse = " + "))), data = d_full)
print(cox.zph(m_full))
pdf("results/TME/schoenfeld_full_model.pdf", width = 8, height = 6)
plot(cox.zph(m_full)[1])
dev.off()
