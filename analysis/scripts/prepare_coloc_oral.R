# Extract the same 22-gene cis windows from the strict oral OLP endpoint, so that
# colocalisation can be evaluated in the third case definition.
# Writes data/processed/coloc/coloc_fg_oral.csv (the existing wide/inclavo files are untouched).
suppressPackageStartupMessages(library(data.table))

eqwin <- fread("data/processed/coloc/coloc_eqtl_windows.csv")
all_rs <- unique(eqwin$snp)
cat("unique window rsids:", length(all_rs), "\n")

path <- "data/raw/finngen_r13/finngen_R13_K11_ORAL_LICHEN_PLANUS.gz"
cat("reading", basename(path), "...\n")
t0 <- Sys.time()
fg <- fread(path, sep = "\t", header = TRUE,
            select = c("rsids", "ref", "alt", "pval", "beta", "sebeta", "af_alt"),
            colClasses = list(character = c("rsids", "ref", "alt")))
cat("  read", nrow(fg), "rows in", round(difftime(Sys.time(), t0, units = "secs"), 1), "s\n")
setnames(fg, "rsids", "snp")
fg <- fg[!is.na(snp) & snp %in% all_rs]
fg[, snp := trimws(snp)]
fg <- fg[!duplicated(snp)]
fg[, outcome := "ORAL"]
cat("window rsids matched:", nrow(fg), "\n")
fwrite(fg, "data/processed/coloc/coloc_fg_oral.csv")
cat("written data/processed/coloc/coloc_fg_oral.csv\n")
