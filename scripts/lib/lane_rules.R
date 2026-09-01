# lane_rules.R
#
# The chr21 significance-lane rule, extracted from
# scripts/04_chr21_lane_assignment.R so it can be unit-tested. This rule
# produces the headline classification of the analysis, and it has already
# broken once in a way the output-inspection pass caught and the test suite
# did not (an unreachable-lanes regression). It lives here now.
#
# The rule is Hunter et al. (2023)'s own classification, applied on the
# ploidy-corrected scale:
#
#   eligible                = not high-repeat, not low-expression, and the
#                             corrected log2FC is estimable
#   passes_magnitude_filter = eligible AND abs(norm_log2FC) >= deviation_lfc
#   sig_lane                = High_repeats / Low_expression first, then
#                             Expected_dosage for eligible genes inside the
#                             deviation threshold, then DE_low / DE_high by
#                             sign for genes that also clear padj < alpha,
#                             then Not_DE_outside_noise.
#
# ORDER MATTERS. Eligibility is tested BEFORE the magnitude filter: high-repeat
# and low-expression genes are excluded because they are too noisy to assess,
# not because they follow expected dosage, so labelling them Expected_dosage
# would be a false claim. Do not reorder the fcase.

#' Assign the chr21 significance lane.
#'
#' @param dt            data.table with columns norm_log2FC, norm_padj,
#'                      low_expr, high_repeat
#' @param alpha         padj threshold for calling a gene DE
#' @param deviation_lfc absolute corrected-log2FC threshold (Hunter's
#'                      FC >= 1.5 cut is log2(1.5))
#' @return the same data.table, modified by reference, with eligible_idx,
#'         passes_magnitude_filter and sig_lane added
assign_sig_lane <- function(dt, alpha, deviation_lfc) {
  stopifnot(data.table::is.data.table(dt),
            all(c("norm_log2FC", "norm_padj", "low_expr", "high_repeat")
                %in% names(dt)),
            is.numeric(alpha), length(alpha) == 1L,
            is.numeric(deviation_lfc), length(deviation_lfc) == 1L)

  dt[, eligible_idx := !low_expr & !high_repeat & !is.na(norm_log2FC)]
  dt[, passes_magnitude_filter :=
       eligible_idx & abs(norm_log2FC) >= deviation_lfc]

  dt[, sig_lane := data.table::fcase(
    high_repeat == TRUE,                                     "High_repeats",
    low_expr == TRUE,                                        "Low_expression",
    passes_magnitude_filter == FALSE,                        "Expected_dosage",
    !is.na(norm_padj) & norm_padj < alpha & norm_log2FC < 0, "DE_low",
    !is.na(norm_padj) & norm_padj < alpha & norm_log2FC > 0, "DE_high",
    default                                                = "Not_DE_outside_noise")]

  dt[]
}
