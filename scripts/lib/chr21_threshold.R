# chr21_threshold.R
#
# Robust outlier test for chr21 ploidy-corrected deviations.
#
# This machinery has no counterpart in Hunter et al. and is a deliberate
# addition, not a departure: their family-of-four design got an effect-size
# filter for free (small n -> large lfcSE -> only large effects cleared
# padj < 0.01). At n = 399 that protection is gone.
#
# The null is estimated from the chr21 genes themselves, not from non-chr21
# genes: ploidy normalization only acts on chr21, so a diploid gene's
# "corrected" log2FC is a different quantity and cannot be the reference.
#
# Scale is MAD rather than SD because the chr21 SD is inflated by the very
# outliers being detected (SD/MAD = 1.49; removing 5 of 160 genes moves SD 19%
# and MAD 1%, and those 5 include OLIG2, a primary candidate).
#
# Callers MUST apply the low-expression and high-repeat filters before
# estimating the null - log2FC variance scales with expression.

MIN_NULL_N <- 10

n_missing <- function(x) sum(is.na(x))

#' Estimate the robust null from chr21 corrected log2 fold changes.
#' @return list(center, scale, n); scale is mad(), already x1.4826.
chr21_null <- function(log2fc) {
  stopifnot(is.numeric(log2fc))
  x <- log2fc[!is.na(log2fc)]
  if (length(x) < MIN_NULL_N) {
    stop("need at least ", MIN_NULL_N, " non-missing values to estimate a null; got ",
         length(x))
  }
  s <- mad(x)
  if (!is.finite(s) || s <= 0) {
    stop("degenerate null scale (mad = ", s, "); the input has no spread")
  }
  list(center = median(x), scale = s, n = length(x))
}

robust_z <- function(log2fc, null) {
  stopifnot(is.list(null), all(c("center", "scale") %in% names(null)))
  if (!is.finite(null$scale) || null$scale <= 0) stop("null scale must be positive")
  (log2fc - null$center) / null$scale
}

#' Two-sided BH q-values. The chr21 bulk is normal once the 10% most extreme
#' are trimmed (Shapiro-Wilk p = 0.17 vs 6e-6 untrimmed), which licenses a
#' normal reference. Re-check that after Task 2 changes the values.
outlier_fdr <- function(z) p.adjust(2 * pnorm(-abs(z)), method = "BH")

effective_k <- function(z, q, alpha) {
  sel <- !is.na(q) & !is.na(z) & q < alpha
  if (!any(sel)) return(NA_real_)
  min(abs(z[sel]))
}

#' Gene counts across k, with the normal-null expectation. Report this instead
#' of a single count: the DE_low set is threshold-driven while DE_high is not.
k_sensitivity <- function(z, k_values = c(1, 1.5, 2, 2.5, 3, 3.5)) {
  z <- z[!is.na(z)]
  n_genes <- vapply(k_values, function(k) sum(abs(z) >= k), integer(1))
  data.frame(k = k_values, n_genes = n_genes,
             pct_genes = 100 * n_genes / length(z),
             expected_pct_normal = 100 * 2 * pnorm(-k_values))
}
