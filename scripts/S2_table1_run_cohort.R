# S2_table1_run_cohort.R
#
# Purpose: Table 1 (cohort characteristics) for a run's cohort AFTER its own
#          exclusions: samples missing a demographic covariate the run
#          adjusts for (Age_at_visit, Sex, BMI, Sample_source), mosaic T21
#          (runs that exclude it), and samples without a cell-type-
#          proportion (CyTOF) match (composition runs). script
#          00_preprocess_data.R's Table 1 describes the shared, pre-run
#          analysis cohort (302 T21 + 95 Control); this describes what is
#          actually left once a specific run's covariates and exclusions are
#          applied - for baseline (no covariates, no exclusions) the two
#          coincide. Reuses scripts/lib/table1.R's row helpers so both
#          tables are built the same way; only the input cohort differs.
#
# Usage:  Rscript scripts/S2_table1_run_cohort.R --run adjusted
# Inputs: results/runs/<name>/processed/cohort_roster.csv (script 01);
#         data/P4C_Comorbidity_020921.tsv (optional)
# Outputs: results/runs/<name>/tables/table1_analysis_cohort.{csv,md}

suppressPackageStartupMessages(library(data.table))
source("scripts/lib/run.R"); run <- load_run()
source("scripts/lib/table1.R")
cat(sprintf("=== Table 1: %s run cohort ===\n\n", run$name))

roster <- fread(run$processed("cohort_roster.csv"))
roster[excluded_reason == "", excluded_reason := NA_character_]
cohort <- roster[is.na(excluded_reason)]
excl   <- roster[!is.na(excluded_reason), .N, by = excluded_reason][order(-N)]

cat(sprintf("  %s run cohort: T21 %d, Control %d, total %d (of %d in the roster)\n",
            run$name, sum(cohort$Karyotype == "T21"), sum(cohort$Karyotype == "Control"),
            nrow(cohort), nrow(roster)))
cat("  Excluded, by reason:\n"); print(as.data.frame(excl))

header <- rbindlist(c(
  list(data.table(characteristic = sprintf("Subjects in %s run cohort", run$name), level = "n",
                   t21 = as.character(sum(cohort$Karyotype == "T21")),
                   control = as.character(sum(cohort$Karyotype == "Control")),
                   overall = as.character(nrow(cohort)),
                   p_value = NA_real_, test = "", n_missing = NA_integer_)),
  lapply(seq_len(nrow(excl)), function(i) {
    ex <- roster[excluded_reason == excl$excluded_reason[i]]
    data.table(characteristic = sprintf("Excluded: %s", excl$excluded_reason[i]), level = "n",
               t21 = as.character(sum(ex$Karyotype == "T21")),
               control = as.character(sum(ex$Karyotype == "Control")),
               overall = as.character(excl$N[i]),
               p_value = NA_real_, test = "", n_missing = NA_integer_)
  })))

tbl1 <- rbindlist(list(header, build_table1_rows(cohort, "data/P4C_Comorbidity_020921.tsv")))
tbl1[, p_value := ifelse(is.na(p_value), "", format.pval(p_value, digits = 2, eps = 1e-4))]

fwrite(tbl1, run$table("table1_analysis_cohort.csv"))
writeLines(c(
  sprintf("# Table 1. Characteristics of the %s run cohort", run$name), "",
  sprintf("Subjects kept in the %s run after its own exclusions (n = %d of %d in the pre-run analysis cohort).",
          run$name, nrow(cohort), nrow(roster)),
  if (length(run$covariates))
    sprintf("Requires every covariate the run adjusts for: %s.", paste(run$covariates, collapse = ", ")),
  if (!is.null(run$composition)) "Requires a cell-type-proportion (CyTOF) match.",
  if (length(run$exclude))
    sprintf("Excludes %s.", paste(sprintf("%s = %s", names(run$exclude), unlist(run$exclude)), collapse = "; ")),
  "Race and ethnicity are not recorded in the available metadata.", "",
  paste("|", paste(names(tbl1), collapse = " | "), "|"),
  paste("|", paste(rep("---", ncol(tbl1)), collapse = " | "), "|"),
  apply(tbl1, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))),
  run$table("table1_analysis_cohort.md"))

for (suffix in c("csv", "md")) {
  f <- run$table(sprintf("table1_analysis_cohort.%s", suffix))
  if (!file.exists(f)) stop("failed to write ", f)
  cat("  Wrote ", f, "\n", sep = "")
}
print(tbl1)
