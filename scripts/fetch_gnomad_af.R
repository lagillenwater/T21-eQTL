# fetch_gnomad_af.R
#
# Purpose: Cache gnomAD v4.1 genomes ALT allele frequencies for the chr21
#          positions of a run's cis-variant universe. Script 02 calls the same
#          helper itself, so this script is only needed to warm the cache
#          ahead of a pipeline run, or to refresh it after the variant
#          universe changes.
#
# The remote sites VCF is 7.8 GB and is never downloaded: the tabix index is
# fetched and only the compressed blocks covering the requested positions are
# transferred (scripts/lib/gnomad.R). Expect a few minutes on a first run and
# no network at all once every position is cached.
#
# Inputs:  results/runs/<name>/processed/eqtl_target_variants.csv (script 02)
#          or, with --positions <file>, one chr21 position per line.
# Outputs: data/gnomad/gnomad_chr21_af.csv
#          data/gnomad/gnomad_chr21_queried_positions.txt
#
# Usage:
#   Rscript scripts/fetch_gnomad_af.R --run baseline
#   Rscript scripts/fetch_gnomad_af.R --positions data/my_positions.txt

suppressPackageStartupMessages(library(data.table))
source("scripts/lib/gnomad.R")

args <- commandArgs(trailingOnly = TRUE)
i    <- match("--positions", args)

if (!is.na(i)) {
  if (length(args) <= i) stop("--positions given without a file", call. = FALSE)
  pos <- as.integer(readLines(args[i + 1]))
  cat(sprintf("=== gnomAD AF cache: %d positions from %s ===\n\n",
              length(unique(pos[!is.na(pos)])), args[i + 1]))
} else {
  source("scripts/lib/run.R"); run <- load_run()
  path <- run$processed("eqtl_target_variants.csv")
  if (!file.exists(path))
    stop("run scripts/02_filter_genotypes.R --run ", run$name, " first - ",
         basename(path), " is missing", call. = FALSE)
  pos <- fread(path)$POS
  cat(sprintf("=== gnomAD AF cache [run: %s]: %d cis positions ===\n\n",
              run$name, length(unique(pos[!is.na(pos)]))))
}

af     <- ensure_gnomad_af(pos)
status <- attr(af, "gnomad_status")

cat(sprintf("\n  status:            %s\n", status))
cat(sprintf("  rows in cache:     %d\n", nrow(read_gnomad_cache()$af)))
cat(sprintf("  positions covered: %d of %d requested\n",
            uniqueN(af$POS), uniqueN(pos[!is.na(pos)])))
cat(sprintf("  cache files:       %s, %s\n", GNOMAD_CACHE_AF, GNOMAD_CACHE_POSITIONS))

if (!identical(status, "ok"))
  stop("gnomAD fetch incomplete: ", status, call. = FALSE)

cat("\n=== gnomAD AF cache complete ===\n")
