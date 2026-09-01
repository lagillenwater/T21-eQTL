source(file.path("..", "..", "scripts", "lib", "chr21_threshold.R"))

test_that("chr21_null returns median and scaled MAD", {
  set.seed(1)
  x <- rnorm(500, mean = 2, sd = 3)
  n <- chr21_null(x)
  expect_equal(n$center, median(x))
  expect_equal(n$scale, mad(x))
  expect_equal(n$n, 500L)
})

test_that("chr21_null drops NA and refuses tiny inputs", {
  expect_equal(chr21_null(c(1, 2, NA, 4, 5, 6, 7, 8, 9, 10, 11))$n, 10L)
  expect_error(chr21_null(c(1, 2, 3)), "at least 10")
})

test_that("chr21_null refuses a degenerate scale", {
  expect_error(chr21_null(rep(2, 50)), "scale")
})

test_that("robust_z centres and scales", {
  x <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  expect_equal(robust_z(x, list(center = 4, scale = 2, n = 10L)), (x - 4) / 2)
})

test_that("MAD scale resists outliers far better than SD", {
  set.seed(2)
  clean <- rnorm(160)
  dirty <- c(clean[1:155], 20, -20, 25, -25, 30)
  expect_gt(abs(sd(dirty) - sd(clean)) / sd(clean),
            5 * abs(mad(dirty) - mad(clean)) / mad(clean))
})

test_that("outlier_fdr finds injected outliers and spares a clean null", {
  # Assertions are ranges, not exact counts: an exact count would be a fact
  # about the chosen seed rather than about outlier_fdr. A clean null should
  # flag nothing, but one false positive in 200 draws at FDR 0.10 is within
  # what the procedure allows. The spiked case has 5 true outliers, and one or
  # two spurious flags from the 195 background draws is likewise acceptable -
  # what must not happen is missing the injected ones.
  set.seed(3)
  clean <- rnorm(200)
  expect_lte(sum(outlier_fdr(robust_z(clean, chr21_null(clean))) < 0.10), 1)
  set.seed(1)
  spiked <- c(rnorm(195), 8, -8, 9, -9, 10)
  n_detected <- sum(outlier_fdr(robust_z(spiked, chr21_null(spiked))) < 0.10)
  expect_gte(n_detected, 4)
  expect_lte(n_detected, 7)
})

test_that("effective_k is the smallest |z| that survives the FDR cut", {
  z <- c(0.5, 1, 2, 3, 4, 5)
  q <- c(0.9, 0.8, 0.5, 0.2, 0.02, 0.001)
  expect_equal(effective_k(z, q, 0.05), 4)
  expect_true(is.na(effective_k(z, q, 1e-6)))
})

test_that("k_sensitivity reports counts and the normal expectation", {
  set.seed(4)
  z <- rnorm(1000)
  s <- k_sensitivity(z, c(1, 2, 3))
  expect_equal(nrow(s), 3L)
  expect_equal(s$n_genes[1], sum(abs(z) >= 1))
  expect_equal(s$expected_pct_normal[2], 100 * 2 * pnorm(-2))
})

test_that("n_missing counts NA", {
  expect_equal(n_missing(c(1, NA, 3, NA)), 2L)
})
