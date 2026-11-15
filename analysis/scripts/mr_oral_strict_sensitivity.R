# Supplementary analysis: FinnGen R13 strict oral OLP endpoint (K11_ORAL_LICHEN_PLANUS,
# 920 cases / 499,266 controls) added as a third, most specific case definition.
#
# The two primary definitions (WIDE, INCLAVO) and their output tables are untouched;
# everything here is written into results/tables/corrected/ as a sensitivity analysis.
#
# Allele alignment is identical to merge_olp_outcomes.R: FinnGen beta is the ALT-allele
# effect, the exposure effect is per the eQTLGen effect allele, so
#   ref == effect_allele ("same")  -> beta_out = -beta
#   alt == effect_allele ("flip")  -> beta_out = +beta
#
# Usage: Rscript mr_oral_strict_sensitivity.R
suppressPackageStartupMessages(library(data.table))

OUT <- "results/tables/corrected"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

top <- fread("data/processed/eqtl/top_cis_eqtl.csv")
top[, snp := trimws(snp)]
ns <- fread("data/processed/eqtl/top_cis_eqtl_nsamples.csv")

path <- "data/raw/finngen_r13/finngen_R13_K11_ORAL_LICHEN_PLANUS.gz"
cat("== ", basename(path), " ==\n", sep = "")
oral <- fread(path, sep = "\t", header = TRUE,
              select = c("rsids", "ref", "alt", "pval", "beta", "sebeta", "af_alt"),
              colClasses = list(character = c("rsids", "ref", "alt")))
cat("  rows read:", nrow(oral), "\n")
setnames(oral, "rsids", "snp")
oral <- oral[!is.na(snp) & snp != ""]
oral[, snp := trimws(snp)]
oral <- oral[!duplicated(snp)]

m <- merge(top, oral, by = "snp", all.x = TRUE)
m[, allele_state := fifelse(ref == effect_allele & alt == other_allele, "same",
                     fifelse(ref == other_allele & alt == effect_allele, "flip", "mismatch"))]
m[, beta_out := fifelse(allele_state == "same", -beta, beta)]
m[, se_out := sebeta]
m[, pval_out := pval]
m[allele_state == "mismatch", c("beta_out", "se_out", "pval_out", "af_alt") := NA]

cat("  instruments matched:", sum(!is.na(m$beta_out)), "/", nrow(m), "\n")
print(table(m$allele_state, useNA = "ifany"))

d <- merge(m, ns, by = c("gene", "snp"), all.x = TRUE)
d <- d[!is.na(beta_out) & !is.na(NrSamples)]
cat("  analysable genes:", uniqueN(d$gene), "\n")

d[, p_eff := fifelse(allele_state == "same", 1 - af_alt, af_alt)]
d[, MAF := pmin(p_eff, 1 - p_eff)]
d[MAF < 0.005, MAF := 0.005]
d[, se_exp := 1 / sqrt(2 * MAF * (1 - MAF) * NrSamples)]
d[, beta_exp := zscore * se_exp]
stopifnot(all(sign(d$beta_exp) == sign(d$zscore)))

d[, beta_MR := beta_out / beta_exp]
d[, se_MR := se_out / abs(beta_exp)]
d[, Z_MR := (beta_out / se_out) * sign(zscore)]
d[, p_MR := 2 * pnorm(-abs(Z_MR))]
d[, OR := exp(beta_MR)]
d[, OR_lo := exp(beta_MR - 1.96 * se_MR)]
d[, OR_hi := exp(beta_MR + 1.96 * se_MR)]

n_tests <- nrow(d)
d[, p_bonf := pmin(p_MR * n_tests, 1)]
d[, p_fdr := p.adjust(p_MR, method = "BH")]
d[, significant_bonf := p_MR < 0.05 / n_tests]
cat("  genes tested:", n_tests, "| Bonferroni:", formatC(0.05 / n_tests, format = "e"), "\n")
cat("  Bonferroni-significant:", sum(d$significant_bonf), "\n")

# standard-error propagation, same formula as the sensitivity analysis
d[, se_prop := sqrt(se_out^2 / beta_exp^2 + beta_out^2 * se_exp^2 / beta_exp^4)]
d[, p_prop := 2 * pnorm(-abs(beta_MR / se_prop))]
d[, bonf_prop := p_prop < 0.05 / n_tests]

cols <- c("gene", "snp", "chr", "pos", "effect_allele", "other_allele", "allele_state",
          "zscore", "NrSamples", "MAF", "beta_exp", "se_exp", "beta_out", "se_out",
          "beta_MR", "se_MR", "Z_MR", "p_MR", "OR", "OR_lo", "OR_hi", "se_prop",
          "p_prop", "p_bonf", "p_fdr", "significant_bonf", "bonf_prop")
setorderv(d, "p_MR")
fwrite(d[, ..cols], file.path(OUT, "mr_olp_oral_results_corrected.csv"))
cat("written", file.path(OUT, "mr_olp_oral_results_corrected.csv"), "\n")

# ---- five candidates across all three definitions -------------------------
cand <- c("GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1")
rows <- list()
for (tag in c("wide", "inclavo")) {
  x <- fread(sprintf("results/tables/corrected/mr_olp_%s_results_corrected.csv", tag))
  x <- x[gene %in% cand]
  rows[[tag]] <- x[, .(cohort = tag, gene, OR, OR_lo, OR_hi, p_MR)]
}
rows[["oral"]] <- d[gene %in% cand, .(cohort = "ORAL_STRICT", gene, OR, OR_lo, OR_hi, p_MR)]
pool <- rbindlist(rows)
wide_fmt <- dcast(pool, gene ~ cohort, value.var = c("OR", "p_MR"))
setcolorder(wide_fmt, c("gene",
                        "OR_wide", "p_MR_wide", "OR_inclavo", "p_MR_inclavo",
                        "OR_ORAL_STRICT", "p_MR_ORAL_STRICT"))
fwrite(wide_fmt, file.path(OUT, "candidates_three_definitions.csv"))
cat("\n--- five candidates across the three definitions ---\n")
print(wide_fmt, digits = 3)

# ---- do the 22 WIDE hits hold in the strict oral definition? --------------
w <- fread("results/tables/corrected/mr_olp_wide_results_corrected.csv")
hits <- w[significant_bonf == TRUE, .(gene, Z_wide = Z_MR, p_wide = p_MR, OR_wide = OR)]
rep <- merge(hits, d[, .(gene, Z_oral = Z_MR, p_oral = p_MR, OR_oral = OR)],
             by = "gene", all.x = TRUE)
rep[, dir_consistent := sign(Z_wide) == sign(Z_oral)]
rep[, nominal := !is.na(p_oral) & dir_consistent & p_oral < 0.05]
setorderv(rep, "p_oral")
fwrite(rep, file.path(OUT, "wide_hits_in_oral_strict.csv"))
cat("\n--- the 22 WIDE hits in the strict oral definition ---\n")
cat("  with oral estimate:", sum(!is.na(rep$p_oral)), "/", nrow(rep), "\n")
cat("  direction-consistent:", sum(rep$dir_consistent, na.rm = TRUE), "\n")
cat("  nominal (dir + P<0.05):", sum(rep$nominal, na.rm = TRUE), "\n")
cat("  Bonferroni in oral:", sum(!is.na(rep$p_oral) & rep$p_oral < 0.05 / n_tests), "\n")
print(rep[!is.na(p_oral)][1:min(8, .N)], digits = 3)
