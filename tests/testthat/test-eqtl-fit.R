source(file.path("..", "..", "scripts", "lib", "eqtl_fit.R"))

test_that("fit_variants matches lm() for a single variant", {
  set.seed(10)
  n <- 200
  g <- rbinom(n, 3, 0.3)
  e <- 0.7 * g + rnorm(n)
  got <- fit_variants(matrix(g, ncol = 1), e)
  ref <- summary(lm(e ~ g))$coefficients["g", ]
  expect_equal(got$slope[1], unname(ref["Estimate"]), tolerance = 1e-8)
  expect_equal(got$se[1],    unname(ref["Std. Error"]), tolerance = 1e-8)
  expect_equal(got$p[1],     unname(ref["Pr(>|t|)"]), tolerance = 1e-8)
})

test_that("fit_variants matches lm() for every column of a matrix", {
  set.seed(11)
  n <- 150
  G <- matrix(rbinom(n * 5, 3, 0.4), ncol = 5)
  e <- 0.4 * G[, 2] - 0.3 * G[, 4] + rnorm(n)
  got <- fit_variants(G, e)
  ref <- vapply(seq_len(5), function(j) unname(coef(lm(e ~ G[, j]))[2]), numeric(1))
  expect_equal(got$slope, ref, tolerance = 1e-8)
  expect_equal(nrow(got), 5L)
})

test_that("fit_variants returns NA for a monomorphic variant", {
  set.seed(12)
  n <- 100
  got <- fit_variants(cbind(rbinom(n, 3, 0.3), rep(2, n)), rnorm(n))
  expect_false(is.na(got$slope[1]))
  expect_true(is.na(got$slope[2]))
  expect_true(is.na(got$p[2]))
})

test_that("perm_min_p is reproducible and one value per permutation", {
  set.seed(13)
  n <- 120
  G <- matrix(rbinom(n * 8, 3, 0.35), ncol = 8)
  e <- rnorm(n)
  a <- perm_min_p(G, e, n_perm = 50, seed = 99)
  expect_equal(a, perm_min_p(G, e, n_perm = 50, seed = 99))
  expect_length(a, 50L)
  expect_true(all(a >= 0 & a <= 1))
})

test_that("perm_min_p null is not concentrated near zero without signal", {
  set.seed(14)
  n <- 200
  G <- matrix(rbinom(n * 3, 3, 0.4), ncol = 3)
  expect_gt(median(perm_min_p(G, rnorm(n), n_perm = 400, seed = 7)), 0.05)
})

test_that("gene_level_p never returns zero and shrinks with strong signal", {
  # obs = 0.5 is beaten by every permutation null (0.1, 0.2, 0.3), so per the
  # documented formula (1 + sum(mp <= obs)) / (length(mp) + 1) this is fully
  # non-significant: (1 + 3) / 4 = 1. (Brief's original expected value of 1/4
  # was an arithmetic error, inconsistent with this same formula validated by
  # the two assertions below.)
  expect_equal(gene_level_p(0.5, c(0.1, 0.2, 0.3)), 1)
  expect_equal(gene_level_p(0.001, rep(0.5, 99)), 1 / 100)
  expect_gt(gene_level_p(0.9, rep(0.5, 99)), 0.9)
})

test_that("a real eQTL beats its own permutation null", {
  set.seed(15)
  n <- 250
  G <- matrix(rbinom(n * 10, 3, 0.3), ncol = 10)
  e <- 0.9 * G[, 3] + rnorm(n)
  obs <- min(fit_variants(G, e)$p, na.rm = TRUE)
  expect_lt(gene_level_p(obs, perm_min_p(G, e, n_perm = 200, seed = 5)), 0.05)
})

test_that("expr_of_gene aligns counts to subject ids and log-transforms", {
  counts <- data.table::data.table(Gene_name = "GENEA", LAB1 = 3, LAB2 = 7)
  meta   <- data.table::data.table(RecordID = c("S1", "S2"), LabID = c("LAB1", "LAB2"))
  expect_equal(expr_of_gene("GENEA", c("S1", "S2"), counts, meta), log2(c(3, 7) + 1))
  expect_true(all(is.na(expr_of_gene("MISSING", "S1", counts, meta))))
})
