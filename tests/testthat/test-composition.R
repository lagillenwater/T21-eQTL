source(file.path("..", "..", "scripts", "lib", "composition.R"))

test_that("partner_shift + composition_verdict: partners shift with gene -> PROGRAM", {
  set.seed(10)
  n_bg <- 200
  gene_names <- paste0("g", seq_len(n_bg))
  lfc_bg <- setNames(rnorm(n_bg, sd = 0.1), gene_names)
  partner_idx <- 1:20
  lfc_bg[partner_idx] <- 2   # these 20 genes shift hard, same direction as gene

  n_ctrl <- 30
  gene_ctrl <- rnorm(n_ctrl)
  L_ctrl <- matrix(rnorm(n_bg * n_ctrl), nrow = n_bg, ncol = n_ctrl,
                    dimnames = list(gene_names, NULL))
  # Make the 20 partner genes near-perfectly correlated with gene_ctrl.
  for (i in partner_idx) L_ctrl[i, ] <- gene_ctrl + rnorm(n_ctrl, sd = 0.01)

  gene_lfc <- 2
  pm <- partner_shift("TARGET", L_ctrl, gene_ctrl, lfc_bg, n_partners = 20)
  expect_equal(pm, 2, tolerance = 1e-6)

  null <- partner_null(lfc_bg, n_partners = 20, n_draw = 2000, seed = 1)
  expect_length(null, 2000)

  p <- mean(sign(gene_lfc) * null >= sign(gene_lfc) * pm)
  verdict <- composition_verdict(gene_lfc, pm, p)
  expect_equal(verdict, "PROGRAM")
})

test_that("composition_verdict: partners at zero -> GENE-SPECIFIC", {
  gene_lfc <- 2
  partner_lfc <- 0
  p <- 0.5   # a zero partner shift is unremarkable against a near-zero null
  expect_equal(composition_verdict(gene_lfc, partner_lfc, p), "GENE-SPECIFIC")
})

test_that("partner_null is reproducible under the same seed and has length n_draw", {
  bg <- rnorm(50)
  x <- partner_null(bg, n_partners = 5, n_draw = 100, seed = 7)
  y <- partner_null(bg, n_partners = 5, n_draw = 100, seed = 7)
  expect_equal(x, y)
  expect_length(x, 100)

  z <- partner_null(bg, n_partners = 5, n_draw = 100, seed = 8)
  expect_false(identical(x, z))
})

test_that("composition_verdict boundary: share exactly 0.5 is PROGRAM, just under is MIXED", {
  expect_equal(composition_verdict(2, 1, 0.01), "PROGRAM")
  expect_equal(composition_verdict(2, 0.99999, 0.01), "MIXED")
})

test_that("composition_verdict: p >= 0.05 is always GENE-SPECIFIC regardless of share", {
  expect_equal(composition_verdict(2, 2, 0.05), "GENE-SPECIFIC")
  expect_equal(composition_verdict(2, 2, 0.2), "GENE-SPECIFIC")
})
