# verify_callset_counts.R
#
# Purpose: Variant-type counts of the run's genotype call sets, as `bcftools
#          stats` makes them (scripts/lib/callset_counts.R), so the Methods
#          counts can be checked against the files the analysis reads
#          (manuscript #13). A record with an SNV allele and an indel allele
#          counts under both, so the two counts overlap.
#
# Usage:   Rscript scripts/verify_callset_counts.R --run adjusted_ctrl_merged
# Inputs:  the run's genotype PASS files (run$genotypes)
# Outputs: results/runs/<name>/tables/callset_counts.csv

suppressPackageStartupMessages(library(data.table))
source("scripts/lib/callset_counts.R")
source("scripts/lib/run.R"); run <- load_run()

cat(sprintf("=== Call-set counts [run: %s] ===\n\n", run$name))
res <- rbindlist(lapply(names(run$genotypes), function(g) {
  d <- fread(run$genotypes[[g]], select = c("REF", "ALT"), showProgress = FALSE)
  cbind(data.table(set = g, file = run$genotypes[[g]]), callset_counts(d$REF, d$ALT))
}))
fwrite(res, run$table("callset_counts.csv"))
print(res)
cat("\nDone.\n")
