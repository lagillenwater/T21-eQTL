# Tests for scripts/lib/allelic_fc.R - allelic fold change at any ploidy and
# its T21-vs-GTEx comparison. The property that matters: the fit returns the
# per-chromosome allelic ratio, so the same aFC is recovered at ploidy 2 and 3.

library(data.table)

sim <- function(afc, ploidy, n = 400, sd = 0.15, seed = 1) {
  set.seed(seed)
  g <- sample(0:ploidy, n, replace = TRUE)
  list(g = g, y = 4 + afc_curve(g, afc, ploidy) + rnorm(n, sd = sd))
}

test_that("afc_curve is zero at dosage 0 and the aFC at full dosage", {
  expect_equal(afc_curve(0, 1.3, 3), 0)
  expect_equal(afc_curve(3, 1.3, 3), 1.3)
  expect_equal(afc_curve(2, -0.7, 2), -0.7)
  expect_equal(afc_curve(1, 1, 2), log2(1.5))          # one REF + one doubled ALT
})

test_that("fit_afc recovers the same aFC at ploidy 2 and ploidy 3", {
  for (a in c(-1.5, -0.3, 0.4, 2.5)) {
    d2 <- sim(a, 2); d3 <- sim(a, 3)
    f2 <- fit_afc(d2$g, d2$y, 2); f3 <- fit_afc(d3$g, d3$y, 3)
    expect_lt(abs(f2$afc - a), 4 * f2$se)
    expect_lt(abs(f3$afc - a), 4 * f3$se)
    expect_false(f3$at_bound)
  }
})

test_that("a straight line on dosage misses a large aFC that fit_afc recovers", {
  # A rare raising allele: most subjects carry 0 or 1 copies, where the curve
  # is steepest, so three times the linear slope overshoots the 0-to-3 span.
  set.seed(2); g <- rbinom(600, 3, 0.2)
  y <- 4 + afc_curve(g, 4, 3) + rnorm(600, sd = 0.05)
  linear <- 3 * unname(coef(lm(y ~ g))[2])
  expect_gt(abs(linear - 4), 0.5)
  expect_lt(abs(fit_afc(g, y, 3)$afc - 4), 0.1)
})

test_that("counting the other allele negates the aFC and keeps its SE", {
  d <- sim(0.8, 3)
  f <- fit_afc(d$g, d$y, 3); r <- fit_afc(3 - d$g, d$y, 3)
  expect_equal(r$afc, -f$afc, tolerance = 1e-5)
  expect_equal(r$se, f$se, tolerance = 1e-3)
})

test_that("fit_afc standard error matches its sampling spread", {
  est <- vapply(1:200, function(s) { d <- sim(0.6, 3, n = 250, seed = s); fit_afc(d$g, d$y, 3)$afc }, numeric(1))
  d <- sim(0.6, 3, n = 250, seed = 1)
  expect_equal(fit_afc(d$g, d$y, 3)$se, sd(est), tolerance = 0.25)
})

test_that("fit_afc returns NA for a monomorphic variant and rejects bad dosage", {
  expect_true(is.na(fit_afc(rep(1, 20), rnorm(20), 3)$afc))
  expect_error(fit_afc(c(0, 1, 2, 4, 1), rnorm(5), 3), "outside")
})

test_that("gtex_afc_at_variant is exact at the lead variant and scales elsewhere", {
  lead <- gtex_afc_at_variant(0.5, 0.05, 0.5, 0.05, 1.2, 0.15)
  expect_equal(lead$gtex_afc, 1.2); expect_equal(lead$gtex_afc_se, 0.15)
  other <- gtex_afc_at_variant(-0.25, 0.10, 0.5, 0.05, 1.2, 0.15)
  expect_equal(other$gtex_afc, -0.6); expect_equal(other$gtex_afc_se, 0.30)
})

test_that("orient_effect flips only where ALT is not the plotted allele", {
  expect_equal(orient_effect(c(0.4, 0.4, 0.4), c(TRUE, FALSE, NA)), c(0.4, -0.4, NA))
})

test_that("effect_difference combines independent standard errors", {
  d <- effect_difference(1.0, 0.3, 0.4, 0.4)
  expect_equal(d$diff, 0.6); expect_equal(d$diff_se, 0.5)
  expect_equal(d$p, 2 * pnorm(-1.2))
})
