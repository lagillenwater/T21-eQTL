# Tests for scripts/lib/af_bound.R - the allele-frequency bound on how far a
# cis-eQTL can move the T21/D21 group ratio. The property that matters: equal
# allele frequencies leave the ratio unmoved whatever the eQTL strength.

library(data.table)

test_that("allelic_r_from_slope inverts the 0-to-ploidy log2 span", {
  expect_equal(allelic_r_from_slope(0), 0)
  expect_equal(allelic_r_from_slope(1 / 3, ploidy = 3), 1)      # doubles at 3 copies
  expect_equal(allelic_r_from_slope(-1 / 3, ploidy = 3), -0.5)
  expect_true(all(allelic_r_from_slope(c(-2, -0.1, 0.1, 2)) > -1))
})

test_that("predicted_group_shift is zero at equal frequencies for any r", {
  r <- c(-0.9, -0.3, 0, 0.5, 4)
  expect_equal(predicted_group_shift(r, 0.3, 0.3), rep(0, 5))
})

test_that("predicted_group_shift matches the closed form and its sign", {
  expect_equal(predicted_group_shift(1, 0.5, 0.25), log2(1.5 / 1.25))
  expect_gt(predicted_group_shift(1, 0.6, 0.4), 0)     # more of a raising allele
  expect_lt(predicted_group_shift(-0.5, 0.6, 0.4), 0)  # more of a lowering allele
})

test_that("predicted_group_shift rejects impossible inputs", {
  expect_error(predicted_group_shift(-1, 0.2, 0.2), "exceed -1")
  expect_error(predicted_group_shift(0.2, 1.2, 0.2), "outside")
})

test_that("the group shift does not depend on which allele is the reference", {
  # Recoding to the other allele (ALT <-> REF, or to minor/major) flips the
  # within-T21 slope and complements both frequencies; the group ratio is a
  # property of the genotypes, so the predicted shift and bounds must not move.
  slope <- c(-0.3, -0.05, 0.12, 0.4); p_case <- c(0.1, 0.45, 0.62, 0.9); p_ref <- c(0.15, 0.4, 0.7, 0.85)
  expect_equal(predicted_group_shift(allelic_r_from_slope(slope), p_case, p_ref),
               predicted_group_shift(allelic_r_from_slope(-slope), 1 - p_case, 1 - p_ref))
  expect_equal(shift_bounds(allelic_r_from_slope(slope), p_ref),
               shift_bounds(allelic_r_from_slope(-slope), 1 - p_ref))
  d <- c(0.5, -0.4, 0.3, -0.2)
  expect_equal(abs(required_case_af(allelic_r_from_slope(slope), p_ref, d) - p_ref),
               abs(required_case_af(allelic_r_from_slope(-slope), 1 - p_ref, d) - (1 - p_ref)))
})

test_that("shift_bounds brackets every achievable case frequency", {
  r <- c(0.8, -0.4); p_ref <- c(0.3, 0.6)
  b <- shift_bounds(r, p_ref)
  grid <- outer(seq(0, 1, by = 0.05), rep(1, 2))
  s <- predicted_group_shift(rep(r, each = nrow(grid)), as.vector(grid),
                             rep(p_ref, each = nrow(grid)))
  s <- matrix(s, ncol = 2)
  expect_equal(b$bound_low,  apply(s, 2, min))
  expect_equal(b$bound_high, apply(s, 2, max))
  expect_true(all(b$bound_low <= 0 & b$bound_high >= 0))
})

test_that("bound_toward picks the end on the deviation's side", {
  expect_equal(bound_toward(c(-0.2, -0.2), c(0.3, 0.3), c(0.5, -0.5)), c(0.3, -0.2))
  expect_equal(bound_toward(0, 0.3, -0.5), 0)
})

test_that("required_case_af inverts predicted_group_shift", {
  r <- c(0.8, -0.4, 0.3); p_ref <- c(0.3, 0.6, 0.5); p_case <- c(0.45, 0.2, 0.7)
  s <- predicted_group_shift(r, p_case, p_ref)
  expect_equal(required_case_af(r, p_ref, s), p_case)
})

test_that("required_case_af is NA when the shift is out of reach", {
  expect_true(is.na(required_case_af(0.2, 0.3, 1)))     # needs AF > 1
  expect_true(is.na(required_case_af(0, 0.3, 0.1)))     # no allelic effect
  expect_true(is.na(required_case_af(0.5, 0.3, -1)))    # needs AF < 0
})

test_that("binomial_af_draws centres on the frequency with binomial spread", {
  dr <- binomial_af_draws(c(0.3, 0.05), n_copies = 822, B = 20000, seed = 7)
  expect_equal(dim(dr), c(2L, 20000L))
  expect_equal(rowMeans(dr), c(0.3, 0.05), tolerance = 0.002)
  expect_equal(apply(dr, 1, sd), sqrt(c(0.3 * 0.7, 0.05 * 0.95) / 822), tolerance = 0.03)
})

test_that("binomial_af_draws returns whole allele counts over the copies", {
  dr <- binomial_af_draws(0.4, n_copies = 1600, B = 50, seed = 1)
  expect_equal(dr * 1600, round(dr * 1600))
  expect_true(all(dr >= 0 & dr <= 1))
})

test_that("binomial_af_draws keeps a fixed allele fixed and an unknown frequency unknown", {
  dr <- binomial_af_draws(c(0, 1, NA), n_copies = c(822, 822, 822), B = 100, seed = 2)
  expect_equal(unique(dr[1, ]), 0)
  expect_equal(unique(dr[2, ]), 1)
  expect_true(all(is.na(dr[3, ])))
})

test_that("binomial_af_draws takes one copy count or one per variant, and is reproducible", {
  a <- binomial_af_draws(c(0.2, 0.6), n_copies = c(822, 1600), B = 10, seed = 3)
  expect_identical(a, binomial_af_draws(c(0.2, 0.6), n_copies = c(822, 1600), B = 10, seed = 3))
  expect_equal(a[2, ] * 1600, round(a[2, ] * 1600))
  expect_error(binomial_af_draws(c(0.2, 0.6), n_copies = c(822, 1600, 10)), "n_copies")
  expect_error(binomial_af_draws(1.2, n_copies = 822), "\\[0, 1\\]")
})

test_that("r2_from_t recovers the regression R-squared from the residual df", {
  set.seed(3); x <- rep(0:3, 25); y <- 0.2 * x + rnorm(100)
  s <- summary(lm(y ~ x))
  expect_equal(r2_from_t(coef(s)[2, 3], 98), s$r.squared)
})

test_that("r2_from_t with covariates gives the partial R-squared of the variant", {
  set.seed(4); z <- rnorm(100); x <- rep(0:3, 25); y <- 0.2 * x + z + rnorm(100)
  t <- summary(lm(y ~ x + z))$coefficients["x", 3]
  partial <- 1 - sum(resid(lm(y ~ x + z))^2) / sum(resid(lm(y ~ z))^2)
  expect_equal(r2_from_t(t, 97), partial)
})
