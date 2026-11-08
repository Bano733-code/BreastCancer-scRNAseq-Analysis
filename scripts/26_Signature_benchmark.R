# ============================================================
# 26_Signature_benchmark.R
# Is the 60-gene signature more than a proliferation score?
#  (1) proliferation component (cycling-state markers) vs MK component,
#      Cox models, C-index, nested likelihood-ratio tests (LRT)
#  (2) comparison with Ki-67 and with the PAM50 proliferation score
#      (mean z-score of the 11 proliferation genes of Parker et al., 2009)
#  (3) 1,000 random 55-gene sets, from all genes and from the 1,000
#      genes most correlated with the proliferation component
# Input : results/TME/gene_signature_FINAL.csv, gene_signature_pathway_based.csv,
#         results/TME/bulk_signature_scores.csv, METABRIC expression file
# Output: results/TME/signature_components.csv, signature_incremental_LRT.csv,
#         signature_random_benchmark.csv, signature_vs_PAM50_prolif.csv
# NOTE  : needs a lot of RAM (the z-scored METABRIC matrix) and several minutes.
# ============================================================
library(survival); library(dplyr); library(data.table)
set.seed(123)

## ---- 0. clinical data (as in script 25) ----
merged <- read.csv("results/TME/bulk_signature_scores.csv") %>% filter(!is.na(OS_MONTHS), !is.na(OS_STATUS))
merged$event_os  <- ifelse(grepl("DECEASED|1:", merged$OS_STATUS), 1, 0)
merged$event_dss <- ifelse(merged$VITAL_STATUS == "Died of Disease", 1, 0)
for (v in c("CLAUDIN_SUBTYPE", "CHEMOTHERAPY", "HORMONE_THERAPY", "RADIO_THERAPY")) merged[[v]] <- factor(merged[[v]])

## ---- 1. METABRIC expression (as in script 24: first row kept for duplicated symbols) ----
expr <- fread("data/metabric/metabric/data_mrna_illumina_microarray.txt", sep = "\t", header = TRUE, data.table = FALSE)
sc_cols  <- setdiff(colnames(expr), c("Hugo_Symbol", "Entrez_Gene_Id"))
expr_mat <- as.matrix(expr[, sc_cols]); rownames(expr_mat) <- expr[[1]]
expr_mat <- expr_mat[!is.na(rownames(expr_mat)) & !duplicated(rownames(expr_mat)), ]
storage.mode(expr_mat) <- "double"
Z <- t(scale(t(expr_mat))); Z <- Z[rowSums(is.na(Z)) == 0, , drop = FALSE]
rm(expr, expr_mat); gc()

sc  <- function(g) colMeans(Z[intersect(g, rownames(Z)), , drop = FALSE])
z1  <- function(x) as.numeric(scale(x))
sig <- read.csv("results/TME/gene_signature_FINAL.csv")$gene
mk  <- read.csv("results/TME/gene_signature_pathway_based.csv")$gene
pro <- setdiff(sig, mk)
cat("signature:", length(sig), "| MK:", length(mk), "| proliferation:", length(pro),
    "\nfound in METABRIC -> all:", sum(sig %in% rownames(Z)), " pro:", sum(pro %in% rownames(Z)), " MK:", sum(mk %in% rownames(Z)), "\n")
sid <- merged$sample_id
cat("score reproduced, r =", cor(sc(sig)[sid], merged$signature_score), "\n")   # 1
merged$S_all <- z1(sc(sig)[sid]); merged$S_pro <- z1(sc(pro)[sid]); merged$S_mk <- z1(sc(mk)[sid])

# PAM50 proliferation genes (Parker et al., 2009): BIRC5, CCNB1, CDC20, CDCA1 (= NUF2), CEP55,
# KNTC2 (= NDC80), MKI67, PTTG1, RRM2, TYMS, UBE2C; either alias is accepted
alts <- list("BIRC5", "CCNB1", "CDC20", c("NUF2", "CDCA1"), "CEP55", c("NDC80", "KNTC2"), "MKI67", "PTTG1", "RRM2", "TYMS", "UBE2C")
pam  <- sapply(alts, function(a) { h <- intersect(a, rownames(Z)); if (length(h)) h[1] else NA })
cat("PAM50 proliferation genes found:", sum(!is.na(pam)), "of 11\n"); pam <- pam[!is.na(pam)]
merged$S_pam  <- z1(sc(pam)[sid])
merged$S_ki67 <- z1(Z["MKI67", sid])
print(round(cor(merged[, c("S_all", "S_pro", "S_mk", "S_pam", "S_ki67")]), 2))

## ---- 2. models ----
base <- c("AGE_AT_DIAGNOSIS", "CLAUDIN_SUBTYPE")
full <- c(base, "LYMPH_NODES_EXAMINED_POSITIVE", "CHEMOTHERAPY", "HORMONE_THERAPY", "RADIO_THERAPY")
eps  <- list(OS = c("OS_MONTHS", "event_os"), DSS = c("OS_MONTHS", "event_dss"))
covs <- list(age_subtype = base, full = full)
fm   <- function(tm, ev, t) as.formula(paste0("Surv(", tm, ",", ev, ") ~ ", paste(t, collapse = " + ")))
cidx <- function(m) unname(concordance(m)$concordance)
lrt  <- function(a, b) { x <- 2 * (b$loglik[2] - a$loglik[2]); k <- length(coef(b)) - length(coef(a))
  data.frame(chisq = round(x, 1), df = k, p = signif(pchisq(x, k, lower.tail = FALSE), 3)) }

mods <- list("Signature (all 55)" = "S_all", "Proliferation genes only" = "S_pro", "MK genes only" = "S_mk",
             "Prolif + MK (both)" = c("S_pro", "S_mk"), "PAM50 proliferation score" = "S_pam", "Ki-67 (MKI67)" = "S_ki67")
comp <- list(); inc <- list(); pamtab <- list()
for (en in names(eps)) for (cn in names(covs)) {
  tm <- eps[[en]][1]; ev <- eps[[en]][2]; cv <- covs[[cn]]
  d  <- merged[complete.cases(merged[, c(tm, ev, cv, "S_all", "S_pro", "S_mk", "S_pam", "S_ki67")]), ]
  m0 <- coxph(fm(tm, ev, cv), data = d)
  F  <- function(t) coxph(fm(tm, ev, c(cv, t)), data = d)
  for (nm in names(mods)) {
    m <- F(mods[[nm]]); ci <- summary(m)$conf.int; t <- mods[[nm]]
    comp[[length(comp) + 1]] <- cbind(
      data.frame(endpoint = en, adjust = cn, N = nrow(d), model = nm,
                 HR_per_SD = paste(sprintf("%s %.2f (%.2f-%.2f)", t, ci[t, "exp(coef)"], ci[t, "lower .95"], ci[t, "upper .95"]), collapse = "; "),
                 Cindex = round(cidx(m), 4), dC_vs_covariates = round(cidx(m) - cidx(m0), 4)),
      setNames(lrt(m0, m), c("LRT_chisq", "LRT_df", "LRT_p")))
  }
  mp <- F("S_pro"); mm <- F("S_mk"); mb <- F(c("S_pro", "S_mk")); ms <- F("S_all"); mpam <- F("S_pam"); mki <- F("S_ki67")
  inc[[length(inc) + 1]] <- cbind(endpoint = en, adjust = cn, rbind(
    cbind(test = "MK genes added to proliferation genes",  lrt(mp, mb)),
    cbind(test = "Proliferation genes added to MK genes",  lrt(mm, mb)),
    cbind(test = "Signature added to PAM50 proliferation", lrt(mpam, F(c("S_pam", "S_all")))),
    cbind(test = "PAM50 proliferation added to signature", lrt(ms,   F(c("S_all", "S_pam")))),
    cbind(test = "Signature added to Ki-67",               lrt(mki,  F(c("S_ki67", "S_all")))),
    cbind(test = "MK genes added to PAM50 proliferation",  lrt(mpam, F(c("S_pam", "S_mk"))))))
  ci <- summary(mpam)$conf.int["S_pam", ]
  pamtab[[length(pamtab) + 1]] <- rbind(
    data.frame(endpoint = en, adjust = cn, test = "PAM50 prolif. alone",
               result = sprintf("HR %.2f (%.2f-%.2f); C %.4f; dC %.4f", ci[1], ci[3], ci[4], cidx(mpam), cidx(mpam) - cidx(m0)),
               chisq = round(2 * (mpam$loglik[2] - m0$loglik[2]), 1), df = 1,
               p = signif(pchisq(2 * (mpam$loglik[2] - m0$loglik[2]), 1, lower.tail = FALSE), 3)),
    data.frame(endpoint = en, adjust = cn, test = "Signature added to PAM50", result = "", lrt(mpam, F(c("S_pam", "S_all")))),
    data.frame(endpoint = en, adjust = cn, test = "PAM50 added to signature",  result = "", lrt(ms,   F(c("S_all", "S_pam")))),
    data.frame(endpoint = en, adjust = cn, test = "MK genes added to PAM50",   result = "", lrt(mpam, F(c("S_pam", "S_mk")))))
}
comp <- do.call(rbind, comp); inc <- do.call(rbind, inc); pamtab <- do.call(rbind, pamtab)
print(comp, row.names = FALSE); print(inc, row.names = FALSE); print(pamtab, row.names = FALSE)
write.csv(comp, "results/TME/signature_components.csv", row.names = FALSE)
write.csv(inc,  "results/TME/signature_incremental_LRT.csv", row.names = FALSE)
write.csv(pamtab, "results/TME/signature_vs_PAM50_prolif.csv", row.names = FALSE)

## ---- 3. 1,000 random 55-gene signatures ----
Sp       <- z1(sc(pro))
pool_all <- setdiff(rownames(Z), sig)
r        <- as.numeric(Z[pool_all, ] %*% Sp) / (ncol(Z) - 1); names(r) <- pool_all
pool_pro <- names(sort(r, decreasing = TRUE))[1:1000]      # genes most correlated with the proliferation component
ng <- sum(sig %in% rownames(Z)); B <- 1000

rb <- function(tm, ev, cv, pool) {
  d  <- merged[complete.cases(merged[, c(tm, ev, cv)]), ]
  m0 <- coxph(fm(tm, ev, cv), data = d)
  one <- function(s) { d$S <- as.numeric(scale(s)); m <- coxph(fm(tm, ev, c(cv, "S")), data = d)
    c(HR = unname(exp(coef(m)["S"])), chisq = 2 * (m$loglik[2] - m0$loglik[2]), C = cidx(m), p = summary(m)$coefficients["S", "Pr(>|z|)"]) }
  obs <- one(sc(sig)[d$sample_id])
  R <- t(replicate(B, one(colMeans(Z[sample(pool, ng), d$sample_id, drop = FALSE]))))
  data.frame(observed_HR = round(obs["HR"], 3), observed_chisq = round(obs["chisq"], 1), observed_C = round(obs["C"], 4),
             random_median_HR = round(median(R[, "HR"]), 3), random_median_chisq = round(median(R[, "chisq"]), 1),
             pct_random_p_lt_0.05 = round(100 * mean(R[, "p"] < 0.05), 1),
             pct_random_beaten_by_signature = round(100 * mean(R[, "chisq"] < obs["chisq"]), 1),
             empirical_p = round((1 + sum(R[, "chisq"] >= obs["chisq"])) / (B + 1), 4))
}
rnd <- list()
for (en in names(eps)) for (pn in c("all_genes", "prolif_correlated"))
  rnd[[length(rnd) + 1]] <- cbind(endpoint = en, null = pn,
                                  rb(eps[[en]][1], eps[[en]][2], base, if (pn == "all_genes") pool_all else pool_pro))
rnd <- do.call(rbind, rnd); print(rnd, row.names = FALSE)
write.csv(rnd, "results/TME/signature_random_benchmark.csv", row.names = FALSE)
