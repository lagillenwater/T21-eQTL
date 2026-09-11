# run.R
#
# Named, configured pipeline runs. A run is a plain R list in
# config/runs/<name>.R (covariates, composition source, cohort exclusions,
# ploidy, classification thresholds). Every script loads it with
#   source("scripts/lib/run.R"); run <- load_run()
# and resolves its outputs under results/runs/<name>/{processed,tables,figures}.
# Shared, run-independent inputs (raw counts, metadata, genotype PASS files,
# GTEx, CyTOF, karyotype subtype files) keep their data/ paths.
#
# The run name comes from `--run <name>` on the command line, then the
# T21_RUN environment variable, then "baseline".

RUN_COVARIATES <- c("Age_at_visit", "Sex", "BMI", "Sample_source")
RUN_THRESHOLDS <- c("alpha_de", "deviation_lfc", "deviation_lfc_t2",
                    "low_expr_basemean", "fdr_gene", "gtex_pval_keep",
                    "outlier_fdr", "alpha_repro")

#' Validate a run config against the fixed schema; returns it or stops.
validate_run_config <- function(cfg) {
  stop_cfg <- function(...) stop("run config: ", sprintf(...), call. = FALSE)
  for (f in c("name", "covariates", "composition", "exclude", "ploidy", "thresholds"))
    if (!f %in% names(cfg)) stop_cfg("missing field %s", f)
  if (!is.character(cfg$name) || length(cfg$name) != 1 ||
      !grepl("^[A-Za-z0-9_-]+$", cfg$name))
    stop_cfg("name must match ^[A-Za-z0-9_-]+$")
  if (!is.character(cfg$covariates)) stop_cfg("covariates must be character")
  unknown <- setdiff(cfg$covariates, RUN_COVARIATES)
  if (length(unknown)) stop_cfg("unknown covariate %s", unknown[1])
  if (!is.null(cfg$composition)) {
    if (!is.list(cfg$composition)) stop_cfg("composition must be a list or NULL")
    for (f in c("source", "path", "reference"))
      if (is.null(cfg$composition[[f]])) stop_cfg("composition missing %s", f)
  }
  if (!is.list(cfg$exclude)) stop_cfg("exclude must be a list")
  for (k in names(cfg$exclude))
    if (!is.character(cfg$exclude[[k]])) stop_cfg("exclude values must be character")
  if (!is.numeric(cfg$ploidy) || length(cfg$ploidy) != 1) stop_cfg("ploidy must be one number")
  if (!is.list(cfg$thresholds)) stop_cfg("thresholds must be a list")
  for (t in RUN_THRESHOLDS)
    if (is.null(cfg$thresholds[[t]])) stop_cfg("missing threshold %s", t)
  extra <- setdiff(names(cfg$thresholds), RUN_THRESHOLDS)
  if (length(extra)) stop_cfg("unknown threshold %s", extra[1])
  cfg
}

#' Output directories for a run.
run_dirs <- function(name, root = ".") {
  base <- file.path(root, "results", "runs", name)
  list(base = base, processed = file.path(base, "processed"),
       tables = file.path(base, "tables"), figures = file.path(base, "figures"))
}

#' Run name from `--run`, then T21_RUN, then "baseline".
parse_run_name <- function(args) {
  i <- match("--run", args)
  if (!is.na(i)) {
    if (length(args) <= i || !nzchar(args[i + 1]) || startsWith(args[i + 1], "--"))
      stop("run config: --run given without a run name", call. = FALSE)
    return(args[i + 1])
  }
  env <- Sys.getenv("T21_RUN", unset = "")
  if (nzchar(env)) return(env)
  "baseline"
}

#' Load a run: parse the name, validate the config, create the directories,
#' and return the config fields plus path and data accessors.
load_run <- function(args = commandArgs(trailingOnly = TRUE), root = ".") {
  name <- parse_run_name(args)
  path <- file.path(root, "config", "runs", paste0(name, ".R"))
  if (!file.exists(path))
    stop("run config: no file for run '", name, "' at ", path, call. = FALSE)
  cfg <- validate_run_config(eval(parse(path, keep.source = FALSE)))
  if (cfg$name != name)
    stop("run config: name '", cfg$name, "' does not match file ", path, call. = FALSE)
  dirs <- run_dirs(name, root)
  for (d in dirs[c("processed", "tables", "figures")])
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

  run <- cfg
  run$root      <- root
  run$dirs      <- dirs
  run$processed <- function(file) file.path(dirs$processed, file)
  run$table     <- function(file) file.path(dirs$tables, file)
  run$figure    <- function(stem) file.path(dirs$figures, stem)
  run$cohort    <- function() data.table::fread(run$processed("analysis_cohort.csv"))
  # The run's expression artifact (script 01): genes x samples, log2-CPM with
  # the run's covariate effects removed and karyotype kept. Row names are
  # gene names; duplicates are allowed and callers take the first match.
  run$expression <- function() {
    dt <- data.table::fread(run$processed("expression_adjusted.csv"))
    M <- as.matrix(dt[, -1]); rownames(M) <- dt[[1]]
    M
  }
  # Non-chr21 genes with mean raw count >= 25 over the given samples (the
  # run cohort): the expressed-gene pool script S1 uses for the cell-fraction
  # reach counts. Restricting to the cohort keeps the pool independent of
  # samples the run does not use.
  run$partner_pool <- function(counts_dt, lab_ids = NULL) {
    ids <- setdiff(names(counts_dt), c("EnsemblID", "Gene_name", "Chr"))
    if (!is.null(lab_ids)) ids <- intersect(ids, lab_ids)
    unname(rowMeans(as.matrix(counts_dt[, ids, with = FALSE])) >= 25 &
             counts_dt$Chr != "chr21")
  }
  run
}
