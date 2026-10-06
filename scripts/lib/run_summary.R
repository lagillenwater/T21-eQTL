# run_summary.R
#
# The run-level statistics behind the negative conclusion (manuscript #34),
# gathered from one run's tables so two runs can be set side by side
# (scripts/10_compare_runs.R):
#   - deviation bound (script 13): how closely T21 allele frequencies match
#     the primary reference, the share of each deviation the observed
#     frequency difference predicts, the share reachable at frequency
#     extremes, and whether any common, within-T21-supported variant could
#     produce the whole deviation within the 99th percentile of the
#     frequency differences the cohort shows;
#   - base rate (script 16): detection among deviating and Expected-dosage
#     genes, Fisher's exact test;
#   - matched effect size (script 18): observed against expected detections
#     at the deviating genes' GTEx effect sizes, the odds ratio, and the
#     median T21/GTEx effect ratio of each set.

#' One run's negative-conclusion statistics.
#'
#' @param bound_genes af_deviation_bound_genes.csv (script 13)
#' @param background  af_reference_background.csv (script 13)
#' @param base_rate   base_rate_summary.csv (script 16)
#' @param detection   detection_vs_effect_size_summary.csv (script 18)
#' @param ratio       top_variant_effect_ratio_summary.csv (script 18)
#' @return data.table(statistic, value, note), one row per statistic
negative_conclusion_summary <- function(bound_genes, background, base_rate, detection, ratio) {
  g   <- bound_genes[!is.na(variant_id)]
  ref <- background[primary == TRUE]
  dev <- base_rate[set == "deviating genes"]
  exd <- base_rate[set == "expected-dosage genes"]
  aud <- detection[model == "GTEx |aFC|"]
  row <- function(statistic, value, note = NA_character_)
    data.table::data.table(statistic = statistic, value = as.numeric(value), note = note)
  rbind(
    row("n_tested_deviating", nrow(g)),
    row("median_delta_af", ref$median_delta, ref$label),
    row("median_abs_delta_af", ref$median_abs_delta, ref$label),
    row("median_abs_delta_af_sampling",
        if ("median_abs_delta_sampling" %in% names(ref)) ref$median_abs_delta_sampling else NA,
        "expected from binomial sampling alone"),
    row("frac_pred_min", min(g$frac_pred, na.rm = TRUE), "observed delta AF, share of the deviation"),
    row("frac_pred_max", max(g$frac_pred, na.rm = TRUE), "observed delta AF, share of the deviation"),
    row("frac_bound_min", min(g$frac_bound, na.rm = TRUE), "T21 AF at 0 or 1, share of the deviation"),
    row("frac_bound_max", max(g$frac_bound, na.rm = TRUE), "T21 AF at 0 or 1, share of the deviation"),
    row("n_frac_bound_ge_1", sum(g$frac_bound >= 1, na.rm = TRUE), "genes reachable only at a frequency extreme"),
    row("n_route_within_q99", sum(g$easiest_req_pctl < 0.99, na.rm = TRUE),
        "genes with a common supported variant needing a delta AF inside the cohort's 99th percentile"),
    row("min_route_percentile", if (all(is.na(g$easiest_req_pctl))) NA else min(g$easiest_req_pctl, na.rm = TRUE)),
    row("r2_within_t21_min", min(g$r2_within_t21, na.rm = TRUE)),
    row("r2_within_t21_max", max(g$r2_within_t21, na.rm = TRUE)),
    row("pct_detected_deviating", dev$pct_detected, sprintf("%d of %d", dev$n_detected, dev$n_tested)),
    row("pct_detected_expected_dosage", exd$pct_detected, sprintf("%d of %d", exd$n_detected, exd$n_tested)),
    row("fisher_p_detection", exd$fisher_p_vs_deviating),
    row("matched_observed", aud$observed),
    row("matched_expected", aud$expected),
    row("matched_p_at_least", aud$p_at_least),
    row("odds_ratio_matched", aud$odds_ratio_deviating, sprintf("95%% CI %.2f-%.2f", aud$or_ci_lo, aud$or_ci_hi)),
    row("median_effect_ratio_deviating", ratio[set == "Deviating", median_ratio]),
    row("median_effect_ratio_expected_dosage", ratio[set == "Expected dosage", median_ratio],
        sprintf("matched on GTEx |aFC|, n = %d", ratio[set == "Expected dosage", n])),
    row("wilcoxon_p_effect_ratio", ratio$wilcoxon_p[1]))
}
