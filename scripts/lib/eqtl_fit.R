# eqtl_fit.R
#
# Vectorized per-variant expression ~ dosage regressions, and the permutation
# machinery for gene-level significance.
#
# Closed-form matrix algebra rather than repeated lm(): the gene-level
# permutation needs n_variants x n_perm fits per gene, far too many for a loop.
# For centred g and e the slope is sum(g*e)/sum(g^2), so all variants and all
# permutations reduce to one matrix product.

center_cols <- function(M) sweep(M, 2, colMeans(M), "-")

#' Fit expression ~ dosage separately for every variant (column) of G.
#' @param G n x m genotype matrix (alt dosage; 0-3 under trisomy).
#' @return data.frame(slope, se, t, p), m rows. Monomorphic variants give NA.
fit_variants <- function(G, e) {
  stopifnot(is.matrix(G), is.numeric(e), nrow(G) == length(e))
  n <- nrow(G)
  if (n < 3) stop("need at least 3 samples")
  Gc  <- center_cols(G)
  ec  <- e - mean(e)
  Sgg <- colSums(Gc^2)
  Sge <- as.vector(crossprod(Gc, ec))
  See <- sum(ec^2)
  slope <- ifelse(Sgg > 0, Sge / Sgg, NA_real_)
  rss   <- See - ifelse(Sgg > 0, Sge^2 / Sgg, 0)
  se    <- ifelse(Sgg > 0, sqrt((rss / (n - 2)) / Sgg), NA_real_)
  tval  <- slope / se
  data.frame(slope = slope, se = se, t = tval, p = 2 * pt(-abs(tval), df = n - 2))
}

#' Minimum p across variants under permutations of expression.
#'
#' Permuting expression breaks the genotype-expression link while preserving
#' genotype LD and the number of variants tested - exactly the multiplicity the
#' "any variant supports the gene" rule ignores.
perm_min_p <- function(G, e, n_perm = 1000, seed = 42) {
  stopifnot(is.matrix(G), nrow(G) == length(e))
  set.seed(seed)
  n <- nrow(G)
  Gc   <- center_cols(G)
  Sgg  <- colSums(Gc^2)
  keep <- Sgg > 0
  if (!any(keep)) return(rep(NA_real_, n_perm))
  Gc  <- Gc[, keep, drop = FALSE]
  Sgg <- Sgg[keep]
  # One n x n_perm matrix of permuted, centred expression; a single matrix
  # product then gives every variant x permutation slope at once.
  E   <- matrix(e[as.vector(replicate(n_perm, sample.int(n)))], nrow = n)
  Ec  <- sweep(E, 2, colMeans(E), "-")
  Sge <- crossprod(Gc, Ec)
  See <- matrix(colSums(Ec^2), nrow = nrow(Sge), ncol = n_perm, byrow = TRUE)
  slope <- Sge / Sgg
  se    <- sqrt(((See - Sge^2 / Sgg) / (n - 2)) / Sgg)
  apply(2 * pt(-abs(slope / se), df = n - 2), 2, min, na.rm = TRUE)
}

#' Gene-level permutation p-value. The +1s keep it strictly positive, which
#' BH-FDR requires.
gene_level_p <- function(min_p_obs, min_p_perm) {
  mp <- min_p_perm[!is.na(min_p_perm)]
  (1 + sum(mp <= min_p_obs)) / (length(mp) + 1)
}

#' Expression vector for one gene, aligned to subject IDs. Counts are keyed by
#' LabID while genotypes are keyed by subject_id, so this maps through metadata.
expr_of_gene <- function(gene_name, subject_ids, counts, meta_t21) {
  row <- counts[Gene_name == gene_name]
  if (nrow(row) == 0) return(rep(NA_real_, length(subject_ids)))
  # Key on subject_id, the join key the genotype tables use. RecordID is a
  # different ID space entirely (INV... vs HTP...) and keying on it returns all
  # NA for every caller in this pipeline.
  subj <- if ("subject_id" %in% names(meta_t21)) {
    as.character(meta_t21$subject_id)
  } else {
    sub("(?<=[0-9])[A-Z][0-9]*$", "", as.character(meta_t21$LabID), perl = TRUE)
  }
  lab_for_subj <- setNames(as.character(meta_t21$LabID), subj)
  labs <- lab_for_subj[as.character(subject_ids)]
  # Index via a plain named vector rather than row[1, j, with = FALSE]: when
  # labs contains NA (unmatched subject_id), data.table's `[.data.table`
  # errors on an NA column index instead of returning NA.
  row_vals <- suppressWarnings(as.numeric(row[1]))
  names(row_vals) <- names(row)
  log2(unname(row_vals[labs]) + 1)
}
