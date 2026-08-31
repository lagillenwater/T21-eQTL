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
