# Tests for scripts/lib/run_summary.R: the run-level statistics behind the
# negative conclusion, gathered from one run's script 13, 16 and 18 tables.

library(data.table)

fixture <- function() list(
  bound_genes = data.table(
    Gene_name = c("G1", "G2", "G3"), variant_id = c("v1", "v2", "v3"),
    frac_pred = c(-0.05, 0.02, 0.04), frac_bound = c(0.2, 1.4, NA),
    easiest_req_pctl = c(NA, 0.995, 1), r2_within_t21 = c(0.01, 0.3, 0.1)),
  background = data.table(
    reference = c("gtex_af", "gnomad_af"), label = c("GTEx whole blood", "gnomAD v4.1 global"),
    n_variants = c(100L, 100L), median_delta = c(-0.001, 0.006),
    median_abs_delta = c(0.012, 0.023), median_abs_delta_sampling = c(0.009, NA), primary = c(TRUE, FALSE)),
  base_rate = data.table(
    set = c("deviating genes", "expected-dosage genes", "expected-dosage genes, baseMean 98-2231 (deviating range)"),
    n_tested = c(10L, 99L, 59L), n_detected = c(8L, 59L, 34L),
    pct_detected = c(80, 59.6, 57.6), fisher_p_vs_deviating = c(NA, 0.311, NA)),
  detection = data.table(
    model = c("GTEx |aFC|", "GTEx |aFC| + baseMean"), observed = c(8L, 8L),
    expected = c(7.4, 7.34), p_at_least = c(0.494, 0.476),
    odds_ratio_deviating = c(1.47, NA), or_ci_lo = c(0.31, NA), or_ci_hi = c(10.53, NA)),
  ratio = data.table(
    set = c("Expected dosage", "Deviating"), n = c(49L, 10L),
    median_ratio = c(0.76, 0.73), wilcoxon_p = c(0.64, 0.64)))

value_of <- function(s, stat) s[statistic == stat, value]

test_that("negative_conclusion_summary reads the primary reference's frequency agreement", {
  s <- do.call(negative_conclusion_summary, fixture())
  expect_equal(value_of(s, "median_abs_delta_af"), 0.012)
  expect_equal(s[statistic == "median_abs_delta_af", note], "GTEx whole blood")
})

test_that("negative_conclusion_summary gives the range of predicted and reachable fractions", {
  s <- do.call(negative_conclusion_summary, fixture())
  expect_equal(value_of(s, "n_tested_deviating"), 3)
  expect_equal(value_of(s, "frac_pred_min"), -0.05)
  expect_equal(value_of(s, "frac_pred_max"), 0.04)
  expect_equal(value_of(s, "frac_bound_max"), 1.4)
  expect_equal(value_of(s, "n_frac_bound_ge_1"), 1)
})

test_that("negative_conclusion_summary counts genes with a route inside the cohort's 99th percentile", {
  f <- fixture()
  expect_equal(value_of(do.call(negative_conclusion_summary, f), "n_route_within_q99"), 0)
  f$bound_genes$easiest_req_pctl[2] <- 0.95
  expect_equal(value_of(do.call(negative_conclusion_summary, f), "n_route_within_q99"), 1)
})

test_that("negative_conclusion_summary carries the base rate, matched detection and effect ratio", {
  s <- do.call(negative_conclusion_summary, fixture())
  expect_equal(value_of(s, "pct_detected_deviating"), 80)
  expect_equal(value_of(s, "pct_detected_expected_dosage"), 59.6)
  expect_equal(value_of(s, "fisher_p_detection"), 0.311)
  expect_equal(value_of(s, "matched_observed"), 8)
  expect_equal(value_of(s, "matched_expected"), 7.4)
  expect_equal(value_of(s, "odds_ratio_matched"), 1.47)
  expect_equal(value_of(s, "median_effect_ratio_deviating"), 0.73)
  expect_equal(value_of(s, "median_effect_ratio_expected_dosage"), 0.76)
})

test_that("negative_conclusion_summary carries the sampling-only expectation, NA when absent", {
  f <- fixture()
  expect_equal(value_of(do.call(negative_conclusion_summary, f), "median_abs_delta_af_sampling"), 0.009)
  f$background[, median_abs_delta_sampling := NULL]
  expect_true(is.na(value_of(do.call(negative_conclusion_summary, f), "median_abs_delta_af_sampling")))
})
