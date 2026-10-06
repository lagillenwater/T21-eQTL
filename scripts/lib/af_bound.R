# af_bound.R
#
# How far can a cis-eQTL move the T21/D21 group ratio? (manuscript review,
# concern 1). The within-T21 eQTL test asks which people with T21 express a
# gene more; a deviation is a property of the group mean. Under an additive
# allelic model each chr21 copy contributes its allele's output, so the group
# ratio is
#
#   T21 / D21 = 1.5 * (1 + r * p_T21) / (1 + r * p_D21)
#
# with r the relative output difference of the ALT allele over REF and p the
# ALT frequency in each group. Equal frequencies leave the ratio at 1.5 however
# strong the eQTL; only a frequency difference moves it. These helpers turn a
# within-T21 slope into r, the observed frequency difference into the shift it
# predicts, and the frequency extremes (0 or 1 in T21) into the largest shift
# the variant could produce at all.

#' Relative ALT-over-REF allelic output from a within-T21 log2 slope.
#'
#' A trisomic subject with g ALT copies expresses (3 - g) + (1 + r) g, so
#' going from 0 to `ploidy` ALT copies multiplies expression by 1 + r. The
#' log-linear fit spans that range as `ploidy * slope`, giving
#' r = 2^(ploidy * slope) - 1. Exact at the dosage endpoints; for the small
#' slopes of the deviating genes the curvature in between is negligible.
#' @param slope log2 expression per ALT copy (within-T21 fit).
#' @param ploidy chr21 copies in the fitted group (3 for T21).
allelic_r_from_slope <- function(slope, ploidy = 3) {
  2^(ploidy * slope) - 1
}

#' Predicted log2 shift of the group ratio from an ALT-frequency difference.
#'
#' log2((1 + r p_case) / (1 + r p_ref)): zero when the frequencies match,
#' whatever r is. Vectorised over every argument.
predicted_group_shift <- function(r, p_case, p_ref) {
  if (any(r <= -1, na.rm = TRUE)) stop("allelic r must exceed -1 (ALT output > 0)")
  if (any(c(p_case, p_ref) < 0 | c(p_case, p_ref) > 1, na.rm = TRUE))
    stop("allele frequency outside [0, 1]")
  log2((1 + r * p_case) / (1 + r * p_ref))
}

#' Largest shifts reachable if the case frequency went to 0 or to 1.
#' @return data.table(bound_low, bound_high): the reachable interval of the
#'   predicted shift, ordered so bound_low <= bound_high.
shift_bounds <- function(r, p_ref) {
  at0 <- predicted_group_shift(r, 0, p_ref)
  at1 <- predicted_group_shift(r, 1, p_ref)
  data.table::data.table(bound_low = pmin(at0, at1), bound_high = pmax(at0, at1))
}

#' Largest reachable shift in the direction a gene deviates.
#'
#' The signed end of the reachable interval that points the same way as the
#' deviation (0 if the interval never leaves zero in that direction).
bound_toward <- function(bound_low, bound_high, deviation) {
  ifelse(deviation >= 0, pmax(bound_high, 0), pmin(bound_low, 0))
}

#' T21 ALT frequency a variant would need to produce a given group shift.
#'
#' Inverts predicted_group_shift for p_case:
#' p_case = ((1 + r p_ref) 2^shift - 1) / r. NA where r is 0 (no effect) or
#' where the required frequency falls outside [0, 1] (the shift is out of the
#' variant's reach at any frequency).
required_case_af <- function(r, p_ref, shift) {
  p <- ((1 + r * p_ref) * 2^shift - 1) / r
  p[!is.finite(p) | p < 0 | p > 1] <- NA_real_
  p
}

#' Binomial draws of ALT frequencies, for intervals.
#'
#' Each frequency is an ALT count over allele copies (3 per T21 participant,
#' 2 per GTEx donor), so its sampling distribution is taken as binomial:
#' each draw is rbinom(n_copies, p) / n_copies. Vectorised as one rbinom call.
#' @param p ALT frequencies (NA allowed: that row stays NA).
#' @param n_copies allele copies behind each frequency: one value, or one per p.
#' @param B draws per frequency.
#' @return length(p) x B matrix of drawn frequencies.
binomial_af_draws <- function(p, n_copies, B = 2000, seed = 1) {
  if (!length(n_copies) %in% c(1L, length(p))) stop("n_copies must be one value or one per frequency")
  if (any(p < 0 | p > 1, na.rm = TRUE)) stop("allele frequency outside [0, 1]")
  n <- rep_len(n_copies, length(p))
  ok <- !is.na(p)
  out <- matrix(NA_real_, length(p), B)
  set.seed(seed)
  if (any(ok)) {
    k <- stats::rbinom(sum(ok) * B, size = rep(n[ok], B), prob = rep(p[ok], B))
    out[ok, ] <- matrix(k / rep(n[ok], B), sum(ok), B)
  }
  out
}

#' Within-group R-squared of a single-variant regression from its t statistic
#' and residual degrees of freedom: n - 2 without covariates, n - 2 - k with
#' k covariates, when it is the variant's partial R-squared.
r2_from_t <- function(t, df) t^2 / (t^2 + df)
