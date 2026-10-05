# allelic_fc.R
#
# Allelic fold change (aFC) on a chromosome of any ploidy, so a T21 effect is
# on the scale GTEx reports (manuscript review, concern 1 follow-up). Used by
# the spike-in power (script 15) and the effect-size comparison (script 18).
#
# aFC is the log2 ratio of the expression of a chromosome carrying allele A to
# one carrying the other allele (Mohammadi et al. 2017; the `afc` column of
# the GTEx eGenes file). Each chromosome copy contributes its own output, so a
# subject with g copies of A out of `ploidy` expresses
#
#   (ploidy - g) + g * 2^aFC        (in units of one other-allele chromosome)
#
# This is where triploidy enters: in T21 one copy of A is one chromosome of
# three, not one of two, so the same aFC gives a smaller log2 step per copy
# and a curve that flattens as g rises. Fitting the curve, rather than a
# straight line on dosage, returns the same quantity GTEx reports for ploidy 2.

AFC_BOUND <- log2(100)

#' Expected log2 expression at dosage g, relative to g = 0.
afc_curve <- function(g, afc, ploidy) {
  log2((ploidy - g) + g * 2^afc) - log2(ploidy)
}

#' Residual sum of squares of the aFC curve at each candidate aFC, with the
#' intercept profiled out, and with `covariates` (if given) profiled out too:
#' the model is then y = intercept + curve(g; aFC) + covariates. Vectorised
#' over `afc`. project_out() is in eqtl_fit.R.
afc_rss <- function(afc, g, y, ploidy, covariates = NULL) {
  H <- log2(outer(ploidy - g, rep(1, length(afc))) + outer(g, 2^afc))
  R <- y - H
  if (is.null(covariates)) return(colSums(sweep(R, 2, colMeans(R))^2))
  colSums(project_out(R, covariates)^2)
}

#' Least-squares aFC of allele A from dosage and log2 expression.
#'
#' Grid search over [-bound, bound] then a local refinement, so a weak or
#' flat association cannot strand the optimiser. The SE is from the curvature
#' of the residual sum of squares at the optimum; NA when the estimate sits on
#' the bound.
#' @param g copies of allele A per subject (0..ploidy).
#' @param y log2 expression.
#' @param covariates optional matrix fitted jointly with the curve (genotype-
#'   side adjustment); rows with a missing covariate are dropped.
#' @return list(afc, se, intercept, n, at_bound)
fit_afc <- function(g, y, ploidy, bound = AFC_BOUND, grid_n = 401L, covariates = NULL) {
  ok <- !is.na(g) & !is.na(y)
  if (!is.null(covariates)) {
    covariates <- as.matrix(covariates)
    stopifnot(nrow(covariates) == length(g))
    ok <- ok & stats::complete.cases(covariates)
    covariates <- covariates[ok, , drop = FALSE]
  }
  g <- g[ok]; y <- y[ok]; n <- length(y)
  k <- if (is.null(covariates)) 0L else covariate_rank(covariates)
  na_out <- list(afc = NA_real_, se = NA_real_, intercept = NA_real_, n = n, at_bound = NA)
  if (n - k < 4 || stats::var(g) == 0) return(na_out)
  if (any(g < 0 | g > ploidy)) stop("dosage outside [0, ", ploidy, "]", call. = FALSE)
  grid <- seq(-bound, bound, length.out = grid_n)
  a0   <- grid[which.min(afc_rss(grid, g, y, ploidy, covariates))]
  step <- grid[2] - grid[1]
  opt  <- stats::optimize(afc_rss, c(max(a0 - step, -bound), min(a0 + step, bound)),
                          g = g, y = y, ploidy = ploidy, covariates = covariates, tol = 1e-8)
  a <- opt$minimum
  at_bound <- abs(a) > bound - 1e-4
  h  <- 1e-3
  r3 <- afc_rss(c(a - h, a, a + h), g, y, ploidy, covariates)
  d2 <- (r3[1] - 2 * r3[2] + r3[3]) / h^2
  se <- if (at_bound || !is.finite(d2) || d2 <= 0) NA_real_ else sqrt(2 * (r3[2] / (n - 2 - k)) / d2)
  list(afc = a, se = se, intercept = mean(y - afc_curve(g, a, ploidy)), n = n, at_bound = at_bound)
}
