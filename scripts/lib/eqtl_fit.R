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
  P <- 2 * pt(-abs(slope / se), df = n - 2)
  P[!is.finite(P)] <- NA_real_
  # A permutation column with no finite p (e.g. a constant expression vector)
  # must yield NA, not the Inf that min(..., na.rm = TRUE) returns.
  apply(P, 2, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
}

#' Gene-level permutation p-value. The +1s keep it strictly positive, which
#' BH-FDR requires.
gene_level_p <- function(min_p_obs, min_p_perm) {
  mp <- min_p_perm[!is.na(min_p_perm)]
  (1 + sum(mp <= min_p_obs)) / (length(mp) + 1)
}

#' Expression for one gene from the run's artifact matrix, aligned to subject ids.
#'
#' The artifact (script 01, processed/expression_adjusted.csv) is already on
#' the log2 scale; nothing is transformed here. Subjects with no matching
#' LabID column get NA.
expr_from_matrix <- function(gene_name, subject_ids, E, meta_t21) {
  if (!gene_name %in% rownames(E)) return(rep(NA_real_, length(subject_ids)))
  lab_for_subj <- setNames(as.character(meta_t21$LabID), as.character(meta_t21$subject_id))
  labs <- lab_for_subj[as.character(subject_ids)]
  row  <- E[match(gene_name, rownames(E)), ]
  out  <- unname(row[labs])
  out[is.na(labs) | !(labs %in% colnames(E))] <- NA_real_
  out
}
