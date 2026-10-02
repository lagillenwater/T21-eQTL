
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

test_that("perm_min_p returns NA, not Inf, for a constant expression vector", {
  set.seed(16)
  G <- matrix(rbinom(400, 3, 0.3), ncol = 4)
  mp <- perm_min_p(G, rep(5, 100), n_perm = 20, seed = 2)
  expect_true(all(is.na(mp)))
  expect_false(any(is.infinite(mp)))
  expect_equal(gene_level_p(0.01, mp), 1)   # empty null -> (1+0)/(0+1)
})

test_that("expr_from_matrix aligns artifact rows to subject ids without re-logging", {
  E <- rbind(GENEA = c(1.5, 2.5, 3.5)); colnames(E) <- c("HTP0001A", "HTP0002B2", "HTP0003A")
  meta <- data.table::data.table(LabID = colnames(E), subject_id = c("HTP0001", "HTP0002", "HTP0003"))
  expect_equal(expr_from_matrix("GENEA", c("HTP0003", "HTP0001", "HTP0999"), E, meta), c(3.5, 1.5, NA))
  expect_true(all(is.na(expr_from_matrix("NOPE", "HTP0001", E, meta))))
})

# --- covariates on both sides (genotype-side adjustment) ----------------------

sim_cov <- function(n, m, seed) {
  set.seed(seed)
  C <- cbind(age = rnorm(n), male = rbinom(n, 1, 0.5), frac = runif(n))
  G <- matrix(rbinom(n * m, 3, 0.3), ncol = m) + round(C[, "frac"] * 0)   # integer dosages
  G[, 1] <- pmin(3, G[, 1] + (C[, "frac"] > 0.6))                         # one variant tracks a covariate
  list(C = C, G = G)
}

test_that("fit_variants with covariates matches lm(e ~ g + C) for every variant", {
  s <- sim_cov(n = 160, m = 4, seed = 21)
  e <- 0.5 * s$G[, 2] + s$C %*% c(0.8, -0.4, 1.5) + rnorm(160)
  got <- fit_variants(s$G, as.vector(e), covariates = s$C)
  ref <- t(vapply(seq_len(4), function(j)
    summary(lm(e ~ s$G[, j] + s$C))$coefficients[2, c(1, 2, 4)], numeric(3)))
  expect_equal(got$slope, unname(ref[, 1]), tolerance = 1e-8)
  expect_equal(got$se,    unname(ref[, 2]), tolerance = 1e-8)
  expect_equal(got$p,     unname(ref[, 3]), tolerance = 1e-8)
})

test_that("fit_variants returns NA for a variant the covariates explain completely", {
  s <- sim_cov(n = 100, m = 2, seed = 22)
  C <- cbind(s$C, s$G[, 2])
  got <- fit_variants(s$G, rnorm(100), covariates = C)
  expect_false(is.na(got$p[1]))
  expect_true(is.na(got$p[2]))
})

test_that("fit_variants refuses missing values when covariates are given", {
  s <- sim_cov(n = 50, m = 2, seed = 23)
  e <- rnorm(50); e[3] <- NA
  expect_error(fit_variants(s$G, e, covariates = s$C), "missing")
})

test_that("min_p_columns is the per-column minimum of fit_variants p", {
  s <- sim_cov(n = 120, m = 5, seed = 24)
  Y <- matrix(rnorm(120 * 3), ncol = 3)
  for (C in list(NULL, s$C)) {
    ref <- apply(Y, 2, function(y) min(fit_variants(s$G, y, covariates = C)$p, na.rm = TRUE))
    expect_equal(min_p_columns(s$G, Y, covariates = C), ref, tolerance = 1e-10)
  }
})

test_that("perm_min_p without covariates is unchanged by the refactor", {
  set.seed(25)
  n <- 90
  G <- matrix(rbinom(n * 6, 3, 0.3), ncol = 6)
  e <- rnorm(n)
  # the pre-refactor computation, inline
  set.seed(31)
  E   <- matrix(e[as.vector(replicate(40, sample.int(n)))], nrow = n)
  Gc  <- center_cols(G); Sgg <- colSums(Gc^2)
  Ec  <- sweep(E, 2, colMeans(E), "-")
  Sge <- crossprod(Gc, Ec)
  See <- matrix(colSums(Ec^2), nrow = nrow(Sge), ncol = 40, byrow = TRUE)
  P   <- 2 * pt(-abs((Sge / Sgg) / sqrt(((See - Sge^2 / Sgg) / (n - 2)) / Sgg)), df = n - 2)
  expect_equal(perm_min_p(G, e, n_perm = 40, seed = 31), apply(P, 2, min))
})

test_that("perm_min_p with covariates permutes the residualised expression", {
  s <- sim_cov(n = 110, m = 4, seed = 26)
  e <- rnorm(110) + s$C[, 1]
  got <- perm_min_p(s$G, e, n_perm = 30, seed = 3, covariates = s$C)
  e_res <- as.vector(project_out(matrix(e), s$C))
  set.seed(3)
  E <- matrix(e_res[as.vector(replicate(30, sample.int(110)))], nrow = 110)
  expect_equal(got, min_p_columns(s$G, E, covariates = s$C), tolerance = 1e-10)
})
