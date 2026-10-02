# Tests for scripts/lib/alleles.R - minor-allele referencing of the cis-eQTL
# stage. The property that matters throughout: re-referencing is a reflection
# of the regressor, so it flips slope signs and leaves every test statistic
# alone.

library(data.table)

# --- alt_is_minor / maf ------------------------------------------------------

test_that("alt_is_minor splits on 0.5 and propagates NA", {
  expect_equal(alt_is_minor(c(0.01, 0.49, 0.5, 0.51, 0.99)),
               c(TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_true(is.na(alt_is_minor(NA_real_)))
})

test_that("alt_is_minor rejects a frequency outside [0, 1]", {
  expect_error(alt_is_minor(c(0.2, 1.4)), "outside \\[0, 1\\]")
  expect_error(alt_is_minor(-0.1), "outside \\[0, 1\\]")
})

test_that("maf_from_af folds the frequency at 0.5", {
  expect_equal(maf_from_af(c(0.1, 0.5, 0.9)), c(0.1, 0.5, 0.1))
})

test_that("af_is_tie flags only a frequency of exactly 0.5", {
  expect_equal(af_is_tie(c(0.5, 0.5 + 1e-6, NA)), c(TRUE, FALSE, FALSE))
})

# --- allele labels -----------------------------------------------------------

test_that("minor_allele_of names the rarer allele and major_allele_of the commoner", {
  ref <- c("C", "A", "G"); alt <- c("G", "T", "C"); af <- c(0.06, 0.72, NA)
  expect_equal(minor_allele_of(ref, alt, af), c("G", "A", NA))
  expect_equal(major_allele_of(ref, alt, af), c("C", "T", NA))
})

test_that("minor_call returns one row per variant with the tie flagged", {
  got <- minor_call(c("C", "A"), c("G", "T"), c(0.06, 0.5))
  expect_equal(nrow(got), 2L)
  expect_equal(got$minor_allele, c("G", "T"))
  expect_equal(got$maf, c(0.06, 0.5))
  expect_equal(got$maf_tie, c(FALSE, TRUE))
})

# --- slope and dosage alignment ---------------------------------------------

test_that("align_slope_to_minor flips only where ALT is the major allele", {
  expect_equal(align_slope_to_minor(c(0.4, 0.4), c(TRUE, FALSE)), c(0.4, -0.4))
})

test_that("align_slope_to_minor returns NA when the minor allele is unknown", {
  expect_true(is.na(align_slope_to_minor(0.4, NA)))
  expect_true(is.na(align_slope_to_minor(NA_real_, TRUE)))
})

test_that("minor_dosage reflects dosage at the ploidy when ALT is major", {
  expect_equal(minor_dosage(0:3, rep(TRUE, 4), 3), as.numeric(0:3))
  expect_equal(minor_dosage(0:3, rep(FALSE, 4), 3), c(3, 2, 1, 0))
  expect_equal(minor_dosage(c(0, 1), c(FALSE, FALSE), 2), c(2, 1))
})

test_that("minor_dosage refuses a dosage outside the ploidy and a bad ploidy", {
  expect_error(minor_dosage(4, TRUE, 3), "outside \\[0, 3\\]")
  expect_error(minor_dosage(-1, TRUE, 3), "outside \\[0, 3\\]")
  expect_error(minor_dosage(1, TRUE, c(2, 3)), "single number")
})

test_that("T21_CHR21_PLOIDY is the trisomic copy number", {
  expect_equal(T21_CHR21_PLOIDY, 3L)
})

# --- the invariance the alignment rests on ----------------------------------

test_that("re-referencing to the minor allele flips the slope and nothing else", {
  set.seed(21)
  n <- 200
  alt <- rbinom(n, 3, 0.75)                       # ALT is the major allele here
  y   <- 0.5 * alt + rnorm(n)
  fit_alt   <- fit_variants(matrix(alt, ncol = 1), y)
  fit_minor <- fit_variants(matrix(minor_dosage(alt, FALSE, 3), ncol = 1), y)
  expect_equal(fit_minor$slope, -fit_alt$slope, tolerance = 1e-10)
  expect_equal(fit_minor$se, fit_alt$se, tolerance = 1e-10)
  expect_equal(fit_minor$p, fit_alt$p, tolerance = 1e-12)
  expect_equal(align_slope_to_minor(fit_alt$slope, FALSE), fit_minor$slope,
               tolerance = 1e-10)
})

test_that("a sign agreement between two slopes survives re-referencing", {
  gtex <- c(0.3, -0.2, 0.4); t21 <- c(0.1, -0.5, -0.2); am <- c(FALSE, FALSE, TRUE)
  expect_equal(sign(align_slope_to_minor(gtex, am)) == sign(align_slope_to_minor(t21, am)),
               sign(gtex) == sign(t21))
})

# --- cohort frequency and agreement -----------------------------------------

test_that("alt_af_from_dosage divides alt copies by the ploidy", {
  expect_equal(alt_af_from_dosage(c(0, 3, 3), 3), 2 / 3)
  expect_equal(alt_af_from_dosage(c(1, NA, 2), 3), 0.5)
  expect_true(is.na(alt_af_from_dosage(c(NA, NA), 3)))
})

test_that("minor_allele_agrees is NA unless both frequencies are known", {
  expect_equal(minor_allele_agrees(c(0.2, 0.6, 0.2, NA), c(0.3, 0.4, 0.8, 0.2)),
               c(TRUE, FALSE, FALSE, NA))
})

# --- recycling --------------------------------------------------------------
# base::ifelse sizes its answer from the test, so a scalar branch inside a
# nested call used to collapse a whole dosage vector to one value and make the
# variant look monomorphic. These pin the recycling down.

test_that("a scalar alt_minor applies to every element rather than collapsing", {
  expect_equal(minor_dosage(c(0, 1, 2, 3), FALSE, 3), c(3, 2, 1, 0))
  expect_equal(minor_dosage(c(0, 1, 2, 3), TRUE, 3), c(0, 1, 2, 3))
  expect_equal(align_slope_to_minor(c(0.1, -0.2, 0.3), FALSE), c(-0.1, 0.2, -0.3))
})

test_that("a scalar frequency labels every REF/ALT pair", {
  expect_equal(minor_allele_of(c("C", "A"), c("G", "T"), 0.8), c("C", "A"))
  expect_equal(major_allele_of(c("C", "A"), c("G", "T"), 0.8), c("G", "T"))
})
