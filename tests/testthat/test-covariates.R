make_comp <- function(n = 60, k = 4, p = 8, seed = 5) {
  set.seed(seed)
  X <- matrix(runif(n * k, 1, 30), n, k, dimnames = list(paste0("S", seq_len(n)), paste0("c", seq_len(k))))
  B <- matrix(rnorm(p * k, sd = 0.2), p, k)
  Y <- B %*% t(X) + matrix(rnorm(p * n, sd = 0.5), p, n)
  rownames(Y) <- paste0("g", seq_len(p)); colnames(Y) <- rownames(X)
  list(X = X, Y = Y)
}

test_that("read_cytof_wide pivots the long table and records unobserved clusters as 0", {
  f <- tempfile(fileext = ".txt")
  writeLines(c("LabID\tSample_type\tCell_cluster_name\tUnits\tValue",
               "HTP0001A\tWBCs\tCD4+ TCM\tpct\t10.5",
               "HTP0001A\tWBCs\tB cells\tpct\t2.0",
               "HTP0002B2\tWBCs\tCD4+ TCM\tpct\t8.0"), f)
  X <- read_cytof_wide(f)
  expect_equal(dim(X), c(2, 2))
  expect_setequal(rownames(X), c("HTP0001A", "HTP0002B2"))
  expect_setequal(colnames(X), c("CD4+ TCM", "B cells"))
  expect_equal(X["HTP0002B2", "B cells"], 0)
  expect_true(is.na(read_cytof_wide(f, na_as_zero = FALSE)["HTP0002B2", "B cells"]))
})

test_that("composition_fractions drops the reference, orders rows, and stays in percent", {
  W <- matrix(c(50, 30, 20, 40, 40, 20), 2, byrow = TRUE,
              dimnames = list(c("s2", "s1"), c("mono", "b", "t")))
  F <- composition_fractions(W, lab_ids = c("s1", "s2"), reference = "mono")
  expect_equal(rownames(F), c("s1", "s2"))
  expect_equal(colnames(F), c("b", "t"))
  expect_equal(unname(F["s1", ]), c(40, 20))
})

test_that("design_formula puts karyotype last and makes fraction names syntactic", {
  f <- design_formula(c("Age_at_visit", "Sex"), c("CD4+ TCM", "naive B"))
  expect_equal(paste(deparse(f), collapse = ""), "~Age_at_visit + Sex + CD4..TCM + naive.B + karyotype")
  expect_equal(deparse(design_formula(character(0), character(0))), "~karyotype")
  # clusters differing only by a sign must not collapse to one term
  f2 <- design_formula(character(0), c("CD27+ B", "CD27- B", "CD8a+ gd T", "CD8a- gd T"))
  expect_equal(length(all.vars(f2)), 5)
  expect_equal(fraction_colnames(c("CD27+ B", "CD27- B")), c("CD27..B", "CD27..B.1"))
})

test_that("adjust_expression removes covariate effects, keeps the karyotype effect, and is identity with no covariates", {
  set.seed(4); n <- 80; g <- rep(c(TRUE, FALSE), each = 40)
  X <- cbind(age = rnorm(n), bmi = rnorm(n))
  Y <- rbind(a = 2 * g + 1.5 * X[, "age"] + rnorm(n, sd = 0.1),
             b = 0 * g - 2 * X[, "bmi"] + rnorm(n, sd = 0.1))
  colnames(Y) <- paste0("s", 1:n)
  A <- adjust_expression(Y, g, X)
  expect_equal(dim(A), dim(Y))
  # exact property: no covariate effect remains once karyotype is in the model
  expect_equal(unname(coef(lm(A["a", ] ~ g + X))[3:4]), c(0, 0), tolerance = 1e-8)
  expect_equal(unname(coef(lm(A["b", ] ~ g + X))[3:4]), c(0, 0), tolerance = 1e-8)
  expect_equal(unname(mean(A["a", g]) - mean(A["a", !g])), 2, tolerance = 0.1)
  expect_equal(adjust_expression(Y, g, NULL), Y)
  expect_equal(adjust_expression(Y, g, X[, 0, drop = FALSE]), Y)
})

test_that("attribution terms sum to the change in the karyotype coefficient", {
  set.seed(8); n <- 100; g <- rep(c(TRUE, FALSE), each = 50)
  X <- cbind(a = rnorm(n) + 2 * g, b = rnorm(n), c = rnorm(n) - g)
  y <- 1 + 0.5 * g + 0.8 * X[, "a"] - 0.3 * X[, "c"] + rnorm(n, sd = 0.2)
  at <- attribution(y, g, X)
  expect_equal(at$term, c("a", "b", "c"))
  un  <- coef(lm(y ~ g))[2]; adj <- coef(lm(y ~ g + X))[2]
  expect_equal(sum(at$contribution), unname(un - adj), tolerance = 1e-10)
  expect_gt(at$contribution[at$term == "a"], 0)
})

test_that("coef_rows matches lm coefficients row by row", {
  d <- make_comp(n = 80, k = 3, p = 6, seed = 11)
  g <- rep(c(TRUE, FALSE), each = 40)
  X <- cbind(1, as.numeric(g), d$X)
  B <- coef_rows(d$Y, X)
  for (i in seq_len(nrow(d$Y))) expect_equal(unname(B[i, ]), unname(coef(lm(d$Y[i, ] ~ g + d$X))))
})
