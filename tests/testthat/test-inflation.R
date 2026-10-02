test_that("p_from_abs_r matches lm() for one regression", {
  set.seed(20)
  n <- 120
  g <- rbinom(n, 3, 0.3)
  e <- 0.2 * g + rnorm(n)
  ref <- summary(lm(e ~ g))$coefficients["g", "Pr(>|t|)"]
  expect_equal(p_from_abs_r(abs(cor(g, e)), n - 2), ref, tolerance = 1e-10)
})

test_that("project_out on both sides gives the covariate-adjusted lm() p", {
  set.seed(21)
  n <- 150
  C <- cbind(rnorm(n), rnorm(n))
  g <- rbinom(n, 3, 0.4) + 0.5 * (C[, 1] > 0)
  e <- 0.3 * g + C %*% c(1, -0.5) + rnorm(n)
  Ys <- standardise_cols(matrix(e), covariates = C)
  Gs <- standardise_cols(matrix(g), covariates = C)
  r  <- abs(crossprod(Ys, Gs))[1, 1]
  ref <- summary(lm(e ~ g + C))$coefficients["g", "Pr(>|t|)"]
  expect_equal(p_from_abs_r(r, n - 2 - ncol(C)), ref, tolerance = 1e-10)
})

test_that("standardise_cols mean-imputes missing values and zeroes constant columns", {
  M <- cbind(c(1, NA, 3, 4), rep(2, 4))
  S <- standardise_cols(M)
  expect_equal(S[2, 1], 0)
  expect_equal(sum(S[, 1]^2), 1)
  expect_true(all(S[, 2] == 0))
})

test_that("histogram lambda and QQ agree with the per-pair computation", {
  set.seed(22)
  n <- 200
  G <- matrix(rbinom(n * 400, 3, 0.3), n)
  Y <- matrix(rnorm(n * 30), n)
  pos <- seq(1e6, 40e6, length.out = 400)
  tss <- c(rep(20e6, 20), rep(NA, 10))      # last 10 genes on another chromosome
  Ys <- standardise_cols(Y); Gs <- standardise_cols(G)
  bins <- log10p_bins(n - 2)
  sc <- distal_scan(Ys, Gs, tss, pos, gene_group = rep(1:2, c(20, 10)),
                    variant_bin = rep(1L, 400), min_dist = 5e6, bins = bins, block = 7L)
  A <- abs(crossprod(Ys, Gs))
  keep <- abs(outer(ifelse(is.na(tss), -Inf, tss), pos, "-")) >= 5e6
  p <- p_from_abs_r(A[keep], n - 2)
  expect_equal(sum(sc$counts), sum(keep))
  expect_equal(sum(sc$counts[, , 2]), 10 * 400)
  expect_equal(hist_lambda(sc$counts, bins), lambda_from_p(p), tolerance = 2e-3)
  # per-gene medians against a direct masked median
  A[!keep] <- NA
  expect_equal(sc$gene_median_r, apply(A, 1, median, na.rm = TRUE))
  expect_equal(sc$gene_n, rowSums(keep))
  qq <- hist_qq(sc$counts, bins)
  expect_equal(qq$observed[1], -log10(min(p)), tolerance = 2e-4)
  expect_equal(qq$expected[1], -log10(0.5 / length(p)))
})

test_that("distal_scan records hits above the threshold only", {
  set.seed(23)
  n <- 100
  G <- matrix(rbinom(n * 50, 3, 0.4), n)
  Y <- cbind(G[, 10] + rnorm(n, sd = 0.3), rnorm(n))
  sc <- distal_scan(standardise_cols(Y), standardise_cols(G), c(NA, NA), seq_len(50),
                    gene_group = c(1L, 1L), variant_bin = rep(1L, 50), min_dist = 0,
                    bins = log10p_bins(n - 2), hit_x = 8)
  expect_true(all(sc$hits$gene == 1L))
  expect_true(10L %in% sc$hits$variant)
})

test_that("lambda_from_abs_r inverts the median p", {
  df <- 270
  r <- 0.05
  expect_equal(lambda_from_abs_r(r, df),
               qchisq(p_from_abs_r(r, df), 1, lower.tail = FALSE) / qchisq(0.5, 1))
  # the null median |r| gives lambda close to 1
  expect_equal(lambda_from_abs_r(qt(0.75, df) / sqrt(df + qt(0.75, df)^2), df), 1, tolerance = 1e-8)
})

test_that("dosage_matrix averages a subject's sample columns and keeps subject order", {
  dt <- data.table::data.table(POS = c(100L, 200L), REF = c("A", "C"), ALT = c("G", "T"),
                               HTP0001A = c(1, 2), HTP0001B2 = c(3, NA),
                               HTP0002A = c(0, 3), HTP0009A = c(1, 1))
  G <- dosage_matrix(dt, c("HTP0002", "HTP0001", "HTP0005"))
  expect_equal(rownames(G), c("HTP0002", "HTP0001", "HTP0005"))
  expect_equal(colnames(G), c("chr21_100_A_G", "chr21_200_C_T"))
  expect_equal(unname(G["HTP0001", ]), c(2, 2))
  expect_equal(unname(G["HTP0002", ]), c(0, 3))
  expect_true(all(is.na(G["HTP0005", ])))
})

test_that("dosage_stats reports call rate and minor allele frequency", {
  G <- cbind(c(0, 3, 3, NA), c(0, 0, 1, 0))
  st <- dosage_stats(G)
  expect_equal(st$call_rate, c(0.75, 1))
  expect_equal(st$alt_af, c(6 / 9, 1 / 12))
  expect_equal(st$maf, c(1 / 3, 1 / 12))
})

test_that("inverse_normal_cols gives a standard normal rank transform per column", {
  M <- cbind(c(5, 1, 3), c(10, 30, 20))
  Z <- inverse_normal_cols(M)
  expect_equal(Z[, 1], qnorm((c(3, 1, 2) - 0.5) / 3))
  expect_equal(Z[, 2], qnorm((c(1, 3, 2) - 0.5) / 3))
})
