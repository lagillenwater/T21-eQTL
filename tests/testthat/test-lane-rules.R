library(data.table)

ALPHA_T <- 0.01
DEV_T   <- log2(1.5)   # 0.5849625

# Minimal lane input: one row per case, all four required columns present.
lane_dt <- function(...) {
  d <- data.table(...)
  if (!"norm_padj"   %in% names(d)) d[, norm_padj := 1e-10]
  if (!"low_expr"    %in% names(d)) d[, low_expr := FALSE]
  if (!"high_repeat" %in% names(d)) d[, high_repeat := FALSE]
  d
}

test_that("a high-repeat gene is never called Expected_dosage, however large its deviation", {
  # High-repeat genes are excluded because their counts are untrustworthy, not
  # because they follow dosage expectation. Calling one Expected_dosage would
  # assert something the data cannot support.
  d <- lane_dt(norm_log2FC = c(0.01, 5, -5), high_repeat = TRUE)
  assign_sig_lane(d, ALPHA_T, DEV_T)
  expect_equal(d$sig_lane, rep("High_repeats", 3))
  expect_false(any(d$eligible_idx))
  expect_false(any(d$passes_magnitude_filter))
})

test_that("a low-expression gene is never called Expected_dosage", {
  d <- lane_dt(norm_log2FC = c(0.01, 5, -5), low_expr = TRUE)
  assign_sig_lane(d, ALPHA_T, DEV_T)
  expect_equal(d$sig_lane, rep("Low_expression", 3))
  expect_false(any(d$eligible_idx))
})

test_that("high_repeat is tested before low_expr, and both before the magnitude filter", {
  d <- lane_dt(norm_log2FC = 5, low_expr = TRUE, high_repeat = TRUE)
  assign_sig_lane(d, ALPHA_T, DEV_T)
  expect_equal(d$sig_lane, "High_repeats")
})

test_that("the magnitude boundary is inclusive: exactly log2(1.5) is DE", {
  d <- lane_dt(norm_log2FC = c(DEV_T, DEV_T - 1e-9, -DEV_T, -DEV_T + 1e-9))
  assign_sig_lane(d, ALPHA_T, DEV_T)
  expect_equal(d$sig_lane, c("DE_high", "Expected_dosage",
                             "DE_low", "Expected_dosage"))
  expect_equal(d$passes_magnitude_filter, c(TRUE, FALSE, TRUE, FALSE))
})

test_that("sign routes a deviating gene to DE_high or DE_low", {
  d <- lane_dt(norm_log2FC = c(1.2, -1.2))
  assign_sig_lane(d, ALPHA_T, DEV_T)
  expect_equal(d$sig_lane, c("DE_high", "DE_low"))
})

test_that("a gene past the magnitude cut but failing padj is Not_DE_outside_noise", {
  d <- lane_dt(norm_log2FC = c(1.2, -1.2, 1.2),
               norm_padj   = c(0.5, 0.011, NA_real_))
  assign_sig_lane(d, ALPHA_T, DEV_T)
  expect_equal(d$sig_lane, rep("Not_DE_outside_noise", 3))
})

test_that("NA norm_log2FC is never DE and never errors", {
  d <- lane_dt(norm_log2FC = c(NA_real_, NA_real_),
               norm_padj   = c(1e-30, NA_real_))
  expect_silent(assign_sig_lane(d, ALPHA_T, DEV_T))
  expect_false(any(d$sig_lane %in% c("DE_high", "DE_low")))
  expect_equal(d$sig_lane, rep("Not_assessable", 2))
  expect_false(any(d$eligible_idx))
  expect_false(any(is.na(d$passes_magnitude_filter)))
})

test_that("every lane in the rule is reachable", {
  # The regression this file exists for was an unreachable lane, caught only by
  # eyeballing the output table. Pin all six.
  d <- lane_dt(
    norm_log2FC = c(1.2, -1.2, 0.1, 1.2, 5,     5),
    norm_padj   = c(1e-9, 1e-9, 1e-9, 0.5, 1e-9, 1e-9),
    low_expr    = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
    high_repeat = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE))
  assign_sig_lane(d, ALPHA_T, DEV_T)
  expect_setequal(d$sig_lane,
                  c("DE_high", "DE_low", "Expected_dosage",
                    "Not_DE_outside_noise", "Low_expression", "High_repeats"))
})

test_that("assign_sig_lane rejects a table missing a required column", {
  d <- data.table(norm_log2FC = 1, norm_padj = 0.001, low_expr = FALSE)
  expect_error(assign_sig_lane(d, ALPHA_T, DEV_T))
})
