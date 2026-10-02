# allelic_fc.R
#
# Allelic fold change (aFC) on a chromosome of any ploidy, and its comparison
# between the T21 cohort and GTEx (manuscript review, concern 1 follow-up: does
# an allele have the same per-copy effect on a trisomic background as in
# euploid blood?).
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
#' intercept profiled out. Vectorised over `afc`.
afc_rss <- function(afc, g, y, ploidy) {
  H <- log2(outer(ploidy - g, rep(1, length(afc))) + outer(g, 2^afc))
  R <- y - H
  colSums(sweep(R, 2, colMeans(R))^2)
}

#' Least-squares aFC of allele A from dosage and log2 expression.
#'
#' Grid search over [-bound, bound] then a local refinement, so a weak or
#' flat association cannot strand the optimiser. The SE is from the curvature
#' of the residual sum of squares at the optimum; NA when the estimate sits on
#' the bound.
#' @param g copies of allele A per subject (0..ploidy).
#' @param y log2 expression.
#' @return list(afc, se, intercept, n, at_bound)
fit_afc <- function(g, y, ploidy, bound = AFC_BOUND, grid_n = 401L) {
  ok <- !is.na(g) & !is.na(y); g <- g[ok]; y <- y[ok]; n <- length(y)
  na_out <- list(afc = NA_real_, se = NA_real_, intercept = NA_real_, n = n, at_bound = NA)
  if (n < 4 || stats::var(g) == 0) return(na_out)
  if (any(g < 0 | g > ploidy)) stop("dosage outside [0, ", ploidy, "]", call. = FALSE)
  grid <- seq(-bound, bound, length.out = grid_n)
  a0   <- grid[which.min(afc_rss(grid, g, y, ploidy))]
  step <- grid[2] - grid[1]
  opt  <- stats::optimize(afc_rss, c(max(a0 - step, -bound), min(a0 + step, bound)),
                          g = g, y = y, ploidy = ploidy, tol = 1e-8)
  a <- opt$minimum
  at_bound <- abs(a) > bound - 1e-4
  h  <- 1e-3
  r3 <- afc_rss(c(a - h, a, a + h), g, y, ploidy)
  d2 <- (r3[1] - 2 * r3[2] + r3[3]) / h^2
  se <- if (at_bound || !is.finite(d2) || d2 <= 0) NA_real_ else sqrt(2 * (r3[2] / (n - 2)) / d2)
  list(afc = a, se = se, intercept = mean(y - afc_curve(g, a, ploidy)), n = n, at_bound = at_bound)
}

#' GTEx aFC at a variant other than the gene's lead variant.
#'
#' GTEx reports aFC only at each gene's lead variant; every other variant has
#' just a slope on normalised expression. Both scales measure the same gene,
#' so the lead variant's own aFC-to-slope ratio converts a slope to aFC units,
#' and its aFC-SE-to-slope-SE ratio converts the SE. Exact at the lead variant
#' itself; an approximation elsewhere (it assumes aFC is proportional to the
#' normalised slope within a gene).
#' @return data.table(gtex_afc, gtex_afc_se)
gtex_afc_at_variant <- function(slope, slope_se, lead_slope, lead_slope_se, lead_afc, lead_afc_se) {
  data.table::data.table(gtex_afc    = slope * lead_afc / lead_slope,
                         gtex_afc_se = slope_se * lead_afc_se / lead_slope_se)
}

#' Re-reference an ALT-coded effect to the allele a panel is drawn on.
#' Swapping the reference allele negates an aFC (and any slope).
orient_effect <- function(effect, alt_is_plot_allele) {
  ifelse(is.na(alt_is_plot_allele), NA_real_, effect * ifelse(alt_is_plot_allele, 1, -1))
}

#' Difference of two independent estimates with a z test.
#' @return data.table(diff, diff_se, diff_lo, diff_hi, z, p)
effect_difference <- function(est1, se1, est2, se2) {
  d  <- est1 - est2
  se <- sqrt(se1^2 + se2^2)
  z  <- d / se
  data.table::data.table(diff = d, diff_se = se, diff_lo = d - stats::qnorm(0.975) * se,
                         diff_hi = d + stats::qnorm(0.975) * se, z = z,
                         p = 2 * stats::pnorm(-abs(z)))
}
