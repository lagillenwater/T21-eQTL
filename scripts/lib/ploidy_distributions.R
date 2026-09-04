# ploidy_distributions.R
#
# Data preparation for scripts/06_chr21_distribution_panel.R: the uncorrected
# and ploidy-corrected log2 fold changes of the same genes, side by side, for
# chromosome 21 (where the correction acts) and a control chromosome (where
# it must not). Extracted so the gene bookkeeping - protein-coding
# restriction, the same gene set on both scales, chromosome naming - is
# unit-tested (tests/testthat/test-ploidy-distributions.R).

#' Long table of log2FC by chromosome and scale.
#'
#' @param res data.frame or data.table from script 01's
#'   deseq2_all_genes_both_analyses.csv: Chr, Gene_type, raw_log2FC,
#'   norm_log2FC. EnsemblID and Gene_name are carried through when present.
#' @param chromosomes Chr values to keep, in display order
#' @param protein_coding_only restrict to Gene_type == "protein_coding"
#' @return data.table(<ids>, chromosome, scale, log2FC). chromosome is a
#'   factor in the order given; scale is a factor: uncorrected,
#'   ploidy-corrected. A gene with an NA estimate on either scale is dropped
#'   from both, so the two curves per chromosome describe the same genes.
#' @section Errors: a requested chromosome that retains no genes stops with
#'   an error, so a Chr naming change ("22" vs "chr22") cannot produce an
#'   empty control curve.
ploidy_distribution_long <- function(res,
                                     chromosomes = c("chr21", "chr22"),
                                     protein_coding_only = TRUE) {
  need <- c("Chr", "Gene_type", "raw_log2FC", "norm_log2FC")
  missing <- setdiff(need, names(res))
  if (length(missing)) {
    stop("res is missing column(s): ", paste(missing, collapse = ", "))
  }
  d <- as.data.table(res)[Chr %in% chromosomes]
  if (protein_coding_only) d <- d[Gene_type == "protein_coding"]
  d <- d[!is.na(raw_log2FC) & !is.na(norm_log2FC)]
  absent <- setdiff(chromosomes, unique(d$Chr))
  if (length(absent)) {
    stop("no genes retained for chromosome(s): ", paste(absent, collapse = ", "))
  }

  id_cols <- intersect(c("EnsemblID", "Gene_name"), names(d))
  long <- melt(d[, c(id_cols, "Chr", "raw_log2FC", "norm_log2FC"), with = FALSE],
               id.vars = c(id_cols, "Chr"),
               measure.vars = c("raw_log2FC", "norm_log2FC"),
               variable.name = "scale", value.name = "log2FC")
  long[, chromosome := factor(Chr, levels = chromosomes)]
  long[, scale := factor(fifelse(scale == "raw_log2FC",
                                 "uncorrected", "ploidy-corrected"),
                         levels = c("uncorrected", "ploidy-corrected"))]
  long[, Chr := NULL]
  setcolorder(long, c(id_cols, "chromosome", "scale", "log2FC"))
  setorder(long, chromosome, scale)
  long[]
}

#' Per chromosome x scale summary: n, median, mean, sd of log2FC.
ploidy_distribution_stats <- function(long) {
  st <- long[, .(n = .N, median = median(log2FC), mean = mean(log2FC),
                 sd = sd(log2FC)),
             by = .(chromosome, scale)]
  setorder(st, chromosome, scale)
  st[]
}

#' Per-gene shift (ploidy-corrected minus uncorrected), summarised per
#' chromosome. On chr21 the correction is a rigid shift of -log2(1.5); on a
#' control chromosome every shift should be 0.
ploidy_shift_stats <- function(long) {
  id_cols <- intersect(c("EnsemblID", "Gene_name"), names(long))
  if (!length(id_cols)) stop("a gene id column is needed to pair the scales")
  lhs  <- paste(c(id_cols, "chromosome"), collapse = " + ")
  wide <- dcast(long, as.formula(paste(lhs, "~ scale")), value.var = "log2FC")
  wide[, shift := `ploidy-corrected` - uncorrected]
  sh <- wide[, .(n = .N, median_shift = median(shift),
                 min_shift = min(shift), max_shift = max(shift)),
             by = chromosome]
  setorder(sh, chromosome)
  sh[]
}
