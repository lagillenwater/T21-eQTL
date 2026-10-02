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

# A column whose sum of squares after centring (and covariate removal) is at
# or below this is treated as constant: monomorphic, or explained entirely by
# the covariates. Any polymorphic dosage column is many orders above it.
SGG_TOL <- 1e-8

#' Remove covariates, with an intercept, from every column of M by least
#' squares (QR, so a rank-deficient covariate set is handled).
#'
#' Projecting the same covariates out of genotype and expression is the
#' Frisch-Waugh-Lovell form of putting them in the model: the slope and its
#' t statistic equal those of e ~ g + covariates, on n - 2 - k residual
#' degrees of freedom. tensorQTL adjusts both sides the same way.
project_out <- function(M, C) qr.resid(qr(cbind(1, C)), as.matrix(M))

#' Number of covariate dimensions beyond the intercept (the k of n - 2 - k).
covariate_rank <- function(C) qr(cbind(1, C))$rank - 1L

# Residualise G and the expression column(s) on the covariates, refusing
# missing values (the projection needs complete subjects), and return the
# residual degrees of freedom of the single-variant fit.
residualise_fit_inputs <- function(G, Y, covariates) {
  n <- nrow(G)
  if (is.null(covariates)) return(list(G = G, Y = Y, df = n - 2L))
  covariates <- as.matrix(covariates)
  stopifnot(nrow(covariates) == n)
  if (anyNA(G) || anyNA(Y) || anyNA(covariates)) {
    stop("missing values with covariates: subset to complete subjects first", call. = FALSE)
  }
  list(G = project_out(G, covariates), Y = project_out(Y, covariates),
       df = n - 2L - covariate_rank(covariates))
}

#' Fit expression ~ dosage separately for every variant (column) of G.
#' @param G n x m genotype matrix (alt dosage; 0-3 under trisomy).
#' @param covariates optional n x k matrix projected out of genotype and
#'   expression alike, so each fit is e ~ g + covariates.
#' @return data.frame(slope, se, t, p), m rows. Monomorphic variants, and
#'   variants the covariates explain completely, give NA.
fit_variants <- function(G, e, covariates = NULL) {
  stopifnot(is.matrix(G), is.numeric(e), nrow(G) == length(e))
  n <- nrow(G)
  if (n < 3) stop("need at least 3 samples")
  r  <- residualise_fit_inputs(G, matrix(e), covariates)
  df <- r$df
  if (df < 1) stop("too few samples for the covariates")
  Gc  <- center_cols(r$G)
  ec  <- as.vector(r$Y) - mean(r$Y)
  Sgg <- colSums(Gc^2)
  ok  <- Sgg > SGG_TOL
  Sge <- as.vector(crossprod(Gc, ec))
  See <- sum(ec^2)
  slope <- ifelse(ok, Sge / Sgg, NA_real_)
  rss   <- See - ifelse(ok, Sge^2 / Sgg, 0)
  se    <- ifelse(ok, sqrt((rss / df) / Sgg), NA_real_)
  tval  <- slope / se
  data.frame(slope = slope, se = se, t = tval, p = 2 * pt(-abs(tval), df = df))
}

#' Smallest p across the variants (columns) of G for every expression column
#' of Y: the vectorised form of min(fit_variants(G, y)$p), one matrix product
#' for all variants and columns. NA for a column with no finite p.
min_p_columns <- function(G, Y, covariates = NULL) {
  stopifnot(is.matrix(G), is.matrix(Y), nrow(G) == nrow(Y))
  r  <- residualise_fit_inputs(G, Y, covariates)
  df <- r$df
  Gc   <- center_cols(r$G)
  Sgg  <- colSums(Gc^2)
  keep <- Sgg > SGG_TOL
  if (!any(keep)) return(rep(NA_real_, ncol(Y)))
  Gc  <- Gc[, keep, drop = FALSE]
  Sgg <- Sgg[keep]
  Yc  <- sweep(r$Y, 2, colMeans(r$Y), "-")
  Sge <- crossprod(Gc, Yc)
  See <- matrix(colSums(Yc^2), nrow = nrow(Sge), ncol = ncol(Y), byrow = TRUE)
  slope <- Sge / Sgg
  se    <- sqrt(((See - Sge^2 / Sgg) / df) / Sgg)
  P <- 2 * pt(-abs(slope / se), df = df)
  P[!is.finite(P)] <- NA_real_
  # A column with no finite p (e.g. a constant expression vector) must yield
  # NA, not the Inf that min(..., na.rm = TRUE) returns.
  apply(P, 2, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
}

#' Minimum p across variants under permutations of expression.
#'
#' Permuting expression breaks the genotype-expression link while preserving
#' genotype LD and the number of variants tested - exactly the multiplicity the
#' "any variant supports the gene" rule ignores. With covariates, the
#' residualised expression is permuted and each permutation gets the same
#' covariate-adjusted fit as the observed data (as FastQTL and tensorQTL do).
perm_min_p <- function(G, e, n_perm = 1000, seed = 42, covariates = NULL) {
  stopifnot(is.matrix(G), nrow(G) == length(e))
  set.seed(seed)
  n <- nrow(G)
  if (!any(colSums(center_cols(G)^2) > SGG_TOL)) return(rep(NA_real_, n_perm))
  if (!is.null(covariates)) e <- as.vector(residualise_fit_inputs(G, matrix(e), covariates)$Y)
  # One n x n_perm matrix of permuted expression; a single matrix product then
  # gives every variant x permutation slope at once.
  E <- matrix(e[as.vector(replicate(n_perm, sample.int(n)))], nrow = n)
  min_p_columns(G, E, covariates)
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
