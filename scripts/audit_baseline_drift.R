# audit_baseline_drift.R
#
# Compare a run's tables against an archived flat results/ tree, table by
# table and cell by cell. Used to prove that the run mechanism reproduces the
# 2026-09-04 outputs before those are archived.
#
# Usage:
#   Rscript scripts/audit_baseline_drift.R --run baseline --archive results/archive/2026-09-04_flat
#
# Pass criteria: every STRICT table identical (numeric tolerance 1e-8) and
# the lane table's classification columns identical. Tables that depend on
# the expression scale (script 03 onward moved from raw counts to the
# log2-CPM artifact) are reported, not gated. Writes <archive>/drift_audit.md
# and exits non-zero on failure.

suppressPackageStartupMessages(library(data.table))
source("scripts/lib/run.R")
args    <- commandArgs(trailingOnly = TRUE)
run     <- load_run(args)
archive <- args[match("--archive", args) + 1]
stopifnot(!is.na(archive), dir.exists(archive))

STRICT <- c("deseq2_all_genes_both_analyses.csv", "deseq2_chr21_genes_both_analyses.csv",
            "deseq2_all_genes_ploidy_normalized.csv", "deseq2_all_genes_no_ploidy_norm.csv",
            "deseq2_chr21_genes_both_analyses_nocooks.csv", "deseq2_cooks_diagnostics.csv",
            "chr21_lane_summary.csv", "chr21_lane_flow.csv",
            "ploidy_correction_distribution_stats.csv", "chr21_k_sensitivity.csv")
LANE_STRICT_COLS <- c("Gene_name", "sig_lane", "tier", "eqtl_lane",
                      "low_expr", "high_repeat")
REPORTED <- c("t21_dosage_per_variant.csv", "eqtl_gene_level_perm.csv",
              "eqtl_negative_controls.csv", "t21_representative_variants.csv")
TOL <- 1e-8
KEYS <- c("Gene_name", "EnsemblID", "variant_id", "sig_lane", "eqtl_lane", "scope",
          "control", "gene_set", "definition", "cluster", "k", "ensembl_stable")

compare_table <- function(name, strict_cols = NULL) {
  fa <- file.path(archive, "tables", name); fb <- run$table(name)
  if (!file.exists(fa)) return(sprintf("- %s: not in archive", name))
  if (!file.exists(fb)) return(sprintf("- %s: MISSING in run", name))
  a <- fread(fa); b <- fread(fb)
  key <- intersect(KEYS, intersect(names(a), names(b)))
  if (length(key)) { setorderv(a, key); setorderv(b, key) }
  # Compare on the columns both tables carry. A strict comparison is limited
  # to strict_cols (a run missing one of them fails); otherwise columns present
  # on one side only are reported, not failed, so a retired annotation column
  # (the 2026-09-10 neighborhood columns) does not mask the row comparison.
  note <- ""
  if (!is.null(strict_cols)) {
    absent <- setdiff(intersect(strict_cols, names(a)), names(b))
    if (length(absent)) return(sprintf("- %s: MISSING column(s) in run: %s", name, paste(absent, collapse = ", ")))
    cols <- intersect(strict_cols, names(a))
  } else {
    only_a <- setdiff(names(a), names(b)); only_b <- setdiff(names(b), names(a))
    if (length(only_a) || length(only_b))
      note <- sprintf(" [columns only in archive: %s; only in run: %s]",
                      if (length(only_a)) paste(only_a, collapse = ", ") else "none",
                      if (length(only_b)) paste(only_b, collapse = ", ") else "none")
    cols <- intersect(names(a), names(b))
  }
  a <- a[, cols, with = FALSE]; b <- b[, cols, with = FALSE]
  if (!identical(dim(a), dim(b)))
    return(sprintf("- %s: DIMENSION %s vs %s%s", name,
                   paste(dim(a), collapse = "x"), paste(dim(b), collapse = "x"), note))
  diffs <- vapply(cols, function(cn) {
    x <- a[[cn]]; y <- b[[cn]]
    both_na <- is.na(x) & is.na(y)
    if (is.numeric(x) && is.numeric(y))
      sum(!both_na & (is.na(x) != is.na(y) | abs(x - y) > TOL), na.rm = TRUE)
    else
      sum(!both_na & (is.na(x) != is.na(y) | as.character(x) != as.character(y)), na.rm = TRUE)
  }, 0)
  if (all(diffs == 0)) sprintf("- %s: identical (%d columns, %d rows)%s", name, length(cols), nrow(a), note)
  else sprintf("- %s: DIFFERS in %s%s", name,
               paste(sprintf("%s (%d cells)", names(diffs)[diffs > 0], diffs[diffs > 0]), collapse = ", "), note)
}

strict_lines <- vapply(STRICT, compare_table, "")
lane_strict  <- compare_table("chr21_lane_assignments.csv", LANE_STRICT_COLS)
lane_all     <- compare_table("chr21_lane_assignments.csv")
reported     <- vapply(REPORTED, compare_table, "")

lines <- c(sprintf("# Drift audit: run `%s` vs `%s` (%s)", run$name, archive, format(Sys.time(), "%Y-%m-%d %H:%M")),
           "", "## Strict tables (must be identical)", strict_lines,
           "", "## Lane table, classification columns (must be identical)", lane_strict,
           "", "## Lane table, all columns (differences listed)", lane_all,
           "", "## Expression-scale-dependent tables (differences listed, expected)", reported)
writeLines(lines, file.path(archive, "drift_audit.md"))
cat(paste(lines, collapse = "\n"), "\n\n")

gated  <- c(strict_lines, lane_strict)
failed <- any(grepl("DIFFERS|DIMENSION|MISSING", gated))
if (failed) { cat("DRIFT AUDIT FAILED\n"); quit(status = 1) } else cat("DRIFT AUDIT PASSED\n")
