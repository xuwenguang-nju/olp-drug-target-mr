# Supplementary coloc.susie run: adds the strict oral case definition and records the
# susie convergence state for every gene x cohort fit.
#
# The canonical results/tables/susie_coloc_results.csv is NOT touched; the output of this
# script goes to results/tables/corrected/. WIDE and INCLAVO are re-fitted here only so that
# niter/converged are captured, and their PP.H4 values are checked against the canonical run
# before the new numbers are used.
#
# Method identical to run_susie_coloc.R: pre-fit each dataset with susie_rss (bounded
# max_iter) and pass the class-"susie" objects to coloc.susie(), priors p1=p2=1e-4, p12=1e-5.
suppressPackageStartupMessages({library(data.table); library(coloc); library(susieR)})

ROOT <- "D:/AI-Wisp-Science/Drug Target MR"
setwd(ROOT)
ld_dir <- file.path(ROOT, "data/processed/coloc/ld_1000g")
OUT <- file.path(ROOT, "results/tables/corrected")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

eqtl <- fread("data/processed/coloc/coloc_eqtl_windows.csv")
fg_list <- list(WIDE = fread("data/processed/coloc/coloc_fg_wide.csv"),
                INCLAVO = fread("data/processed/coloc/coloc_fg_inclavo.csv"),
                ORAL = fread("data/processed/coloc/coloc_fg_oral.csv"))
genes <- unique(eqtl$gene)

meta <- list(WIDE = list(s = 9526 / (9526 + 490660), N = 9526 + 490660),
             INCLAVO = list(s = 2491 / (2491 + 497695), N = 2491 + 497695),
             ORAL = list(s = 920 / (920 + 499266), N = 920 + 499266))
LMAX <- 15
MAX_ITER <- 1000
PRIORS <- list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)

align_sign <- function(effect_allele, ref, alt) {
  fifelse(effect_allele == alt, 1, fifelse(effect_allele == ref, -1, 0))
}

read_ld <- function(g) {
  snps <- fread(file.path(ld_dir, sprintf("%s_snps.csv", g)))
  ld_dt <- fread(file.path(ld_dir, sprintf("%s_ld.csv.gz", g)))
  rs <- ld_dt[[1]]
  LD <- as.matrix(ld_dt[, -1])
  rownames(LD) <- rs; colnames(LD) <- rs
  LD <- (LD + t(LD)) / 2
  list(snps = snps, LD = LD)
}

fit_susie <- function(d, label) {
  z <- d$beta / sqrt(d$varbeta)
  names(z) <- d$snp
  LDsub <- d$LD[d$snp, d$snp]
  fit <- tryCatch(susieR::susie_rss(z = z, R = LDsub, n = d$N, L = LMAX, max_iter = MAX_ITER),
                  error = function(e) e)
  if (inherits(fit, "error")) return(list(err = conditionMessage(fit)))
  fit <- coloc:::annotate_susie(fit, d$snp, LDsub)
  cs <- fit$sets$cs
  n <- if (is.null(cs)) 0 else length(cs)
  list(fit = fit, n_cs = n, niter = fit$niter, converged = fit$converged, err = "")
}

run_cohort <- function(g, ld, fg, meta_one, cohort) {
  ld_snps <- ld$snps; LD <- ld$LD
  eq_g <- eqtl[gene == g][, .(snp, effect_allele, other_allele, zscore, pvalue, NrSamples)]
  fg_g <- fg[, .(snp, ref, alt, beta, sebeta, pval)]
  m <- merge(ld_snps, eq_g, by.x = "rsid", by.y = "snp")
  m <- merge(m, fg_g, by.x = "rsid", by.y = "snp")
  m <- m[!is.na(zscore) & !is.na(beta) & !is.na(sebeta) & pval < 1]
  m[, sign_e := align_sign(effect_allele, ref.x, alt.x)]
  m[, sign_f := align_sign(alt.y, ref.x, alt.x)]
  m <- m[sign_e != 0 & sign_f != 0]
  setorder(m, rsid)
  if (nrow(m) < 10) {
    return(data.table(gene = g, cohort = cohort, n_snps = nrow(m), err = "nsnp<10"))
  }
  N_eqtl <- round(median(m$NrSamples))
  d1 <- list(snp = m$rsid, beta = m$sign_e * m$zscore, varbeta = rep(1, nrow(m)),
             MAF = m$maf, N = N_eqtl, type = "quant", LD = LD)
  d2 <- list(snp = m$rsid, beta = m$sign_f * m$beta, varbeta = m$sebeta^2,
             type = "cc", s = meta_one$s, N = meta_one$N, LD = LD)
  s1 <- fit_susie(d1, "eqtl")
  s2 <- fit_susie(d2, cohort)
  base <- data.table(gene = g, cohort = cohort, n_snps = nrow(m),
                     n_cs_eqtl = s1$n_cs, n_cs_fg = s2$n_cs,
                     niter_eqtl = s1$niter, conv_eqtl = s1$converged,
                     niter_fg = s2$niter, conv_fg = s2$converged)
  if (s1$err != "" || s2$err != "") {
    return(cbind(base, data.table(err = paste(s1$err, s2$err, sep = "; "))))
  }
  res <- tryCatch(do.call(coloc.susie, c(list(s1$fit, s2$fit), PRIORS)),
                  error = function(e) e)
  if (inherits(res, "error")) {
    return(cbind(base, data.table(err = paste0("coloc.susie: ", conditionMessage(res)))))
  }
  if (is.data.table(res) && all(is.na(res$nsnps))) {
    return(cbind(base, data.table(err = "no CS in one/both datasets")))
  }
  summ <- res$summary
  best <- summ[which.max(PP.H4.abf)]
  cbind(base, data.table(n_pairs = nrow(summ),
                         pp_h4_max = round(best$PP.H4.abf, 4),
                         pp_h4_2nd = round(sort(summ$PP.H4.abf, decreasing = TRUE)[2], 4),
                         hit1 = best$hit1, hit2 = best$hit2,
                         err = ""))
}

rows <- list()
for (g in genes) {
  ld <- tryCatch(read_ld(g), error = function(e) e)
  if (inherits(ld, "error")) {
    cat("LD missing for", g, "\n"); next
  }
  for (cf in names(fg_list)) {
    cat(sprintf("  %-9s %-8s ... ", g, cf)); flush.console()
    r <- tryCatch(run_cohort(g, ld, fg_list[[cf]], meta[[cf]], cf),
                  error = function(e) data.table(gene = g, cohort = cf,
                                                 err = conditionMessage(e)))
    cat(if (is.null(r$err) || r$err[1] == "") sprintf("PP.H4 %s\n", r$pp_h4_max[1]) else
        sprintf("ERR %s\n", substr(r$err[1], 1, 50)))
    rows[[paste(g, cf)]] <- r
  }
}

out <- rbindlist(rows, fill = TRUE)
fwrite(out, file.path(OUT, "susie_coloc_three_definitions.csv"))
cat("\nwritten", file.path(OUT, "susie_coloc_three_definitions.csv"), "\n")
cat("fits:", nrow(out), "| converged (fg):", sum(out$conv_fg %in% TRUE, na.rm = TRUE),
    "| converged (eqtl):", sum(out$conv_eqtl %in% TRUE, na.rm = TRUE), "\n")
print(out[gene %in% c("GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1"),
          .(gene, cohort, n_cs_eqtl, n_cs_fg, pp_h4_max, err)], digits = 3)
