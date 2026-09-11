# cohort.R
#
# Analysis-cohort definition.
#
# The cohort is deliberately ASYMMETRIC. Genotypes are used only for the
# within-T21 dosage regressions, which controls never enter, so requiring WGS of
# a control would discard 89 of 95 (94%) for no analytic gain:
#
#   T21     - needs RNA-seq AND WGS  -> 302 of 304
#   Control - needs RNA-seq only     -> 95 of 95
#
# Dropping the 2 ungenotyped T21 makes the DE cohort and the eQTL cohort the
# same people. Previously DE ran on 304 and eQTL on 302.

VCF_FIXED_COLS <- c("CHROM", "POS", "ID", "REF", "ALT",
                    "QUAL", "FILTER", "INFO", "FORMAT")

#' Strip the trailing visit suffix from a LabID to get a subject ID.
#' Matches the rule in scripts/02_filter_genotypes.R (~line 186), with one
#' refinement: the visit letter must follow a digit (lookbehind), so a bare
#' LabID with no visit suffix (e.g. "HTP0003") is left untouched rather than
#' having its final digit-run-preceded letter ("P0003") stripped. Real HTP
#' LabIDs always carry a visit suffix after the numeric subject ID, so this
#' produces identical results to script 02 on real data.
#' Joining genotypes on RecordID instead gives zero matches.
subject_id_from_labid <- function(labid) sub("(?<=[0-9])[A-Z][0-9]*$", "", labid, perl = TRUE)

#' Subject IDs carrying chr21 genotypes, from VCF-style CSV HEADERS ONLY.
#'
#' Reads with nrows = 0 so this never touches the 6 GB body - it must stay cheap
#' enough to run in script 00, before any genotype processing.
wgs_subjects <- function(vcf_paths) {
  cols <- unlist(lapply(vcf_paths, function(f) {
    if (!file.exists(f)) {
      warning("genotype file not found, skipping: ", f)
      return(character(0))
    }
    names(data.table::fread(f, nrows = 0))
  }), use.names = FALSE)
  unique(subject_id_from_labid(setdiff(cols, VCF_FIXED_COLS)))
}

#' Subset metadata to the analysis cohort: genotyped T21 plus ALL controls.
#'
#' Written defensively so it works on either a data.frame/tibble or a
#' data.table: base-R `[` row subsetting with a logical vector works
#' identically for both, unlike data.table's NSE `[expr]` form which only
#' works on a data.table.
analysis_cohort <- function(meta) {
  if (!"has_wgs" %in% names(meta)) {
    stop("metadata needs a has_wgs column; call wgs_subjects() first")
  }
  keep <- (meta$Karyotype == "T21" & meta$has_wgs) | meta$Karyotype != "T21"
  meta[keep, ]
}

#' Run-specific cohort with an excluded_reason per row (NA = kept).
#'
#' Reasons are applied in order and the first failing one is recorded:
#' keep rule (genotyped T21 plus all controls), run exclusions
#' (`run$exclude`, key = metadata column, values to drop), missing run
#' covariates, and, when the run has a composition source, absence from the
#' CyTOF table. The full roster is returned so exclusions stay visible.
build_cohort <- function(meta, run, cytof_ids = NULL) {
  m <- data.table::as.data.table(meta)
  if (!"has_wgs" %in% names(m)) stop("metadata needs a has_wgs column; call wgs_subjects() first")
  reason <- rep(NA_character_, nrow(m))
  keep <- (m$Karyotype == "T21" & m$has_wgs) | m$Karyotype != "T21"
  reason[!keep] <- "keep_rule:no_wgs"
  for (k in names(run$exclude)) {
    if (!k %in% names(m)) stop("exclusion column not in metadata: ", k)
    hit <- is.na(reason) & m[[k]] %in% run$exclude[[k]]
    reason[hit] <- sprintf("exclude:%s=%s", k, m[[k]][hit])
  }
  for (cv in run$covariates) {
    if (!cv %in% names(m)) stop("covariate column not in metadata: ", cv)
    hit <- is.na(reason) & is.na(m[[cv]])
    reason[hit] <- paste0("missing:", cv)
  }
  if (!is.null(run$composition)) {
    hit <- is.na(reason) & !(m$LabID %in% cytof_ids)
    reason[hit] <- "missing:cytof"
  }
  m[, excluded_reason := reason]
  m
}
