
# Builds a partner pool in which the first `n_block` genes form a tight
# co-expression module with the (external) target gene and all shift by
# `block_lfc`; everything else is noise with near-zero log2FC.
make_pool <- function(n_bg, n_block, n_ctrl, block_lfc, seed) {
  set.seed(seed)
  gene_names <- paste0("g", seq_len(n_bg))
  lfc_bg <- setNames(rnorm(n_bg, sd = 0.1), gene_names)
  lfc_bg[seq_len(n_block)] <- block_lfc

  gene_ctrl <- rnorm(n_ctrl)
  L_ctrl <- matrix(rnorm(n_bg * n_ctrl), nrow = n_bg, ncol = n_ctrl,
                   dimnames = list(gene_names, NULL))
  # The module genes are near-perfectly correlated with the target gene, and
  # therefore with each other - which is exactly what the matched null has to
  # reproduce when it seeds inside the module.
  block <- seq_len(n_block)
  L_ctrl[block, ] <- matrix(rep(gene_ctrl, each = n_block),
                            nrow = n_block, ncol = n_ctrl) +
                     matrix(rnorm(n_block * n_ctrl, sd = 0.01),
                            nrow = n_block, ncol = n_ctrl)
  list(L_ctrl = L_ctrl, lfc = lfc_bg, gene_ctrl = gene_ctrl)
}

test_that("partner_shift + composition_verdict: partners shift with gene -> PROGRAM", {
  # 20 partners inside a 1000-gene pool: a matched null draw only reproduces
  # the observed shift when its seed lands in the module (2% of draws), so the
  # module shift stays significant against the harder, matched null.
  d <- make_pool(n_bg = 1000, n_block = 20, n_ctrl = 30,
                 block_lfc = 2, seed = 10)
  gene_lfc <- 2

  pm <- partner_shift("TARGET", d$L_ctrl, d$gene_ctrl, d$lfc, n_partners = 20)
  expect_equal(pm, 2, tolerance = 1e-6)

  null <- partner_null(d$L_ctrl, d$lfc, n_partners = 20, n_draw = 300, seed = 1)
  expect_length(null, 300)

  p <- partner_p(gene_lfc, pm, null)
  expect_lt(p, 0.05)
  expect_equal(composition_verdict(gene_lfc, pm, p), "PROGRAM")
})

test_that("composition_verdict: partners at zero -> GENE-SPECIFIC", {
  # A gene that moves on its own: its partners sit at the background, so the
  # matched null reproduces their shift constantly and p is unremarkable.
  d <- make_pool(n_bg = 1000, n_block = 20, n_ctrl = 30,
                 block_lfc = 0, seed = 11)
  gene_lfc <- 2

  pm <- partner_shift("TARGET", d$L_ctrl, d$gene_ctrl, d$lfc, n_partners = 20)
  expect_lt(abs(pm), 0.1)

  null <- partner_null(d$L_ctrl, d$lfc, n_partners = 20, n_draw = 300, seed = 1)
  p <- partner_p(gene_lfc, pm, null)
  expect_gt(p, 0.05)
  expect_equal(composition_verdict(gene_lfc, pm, p), "GENE-SPECIFIC")
})

test_that("the matched null is far wider than an independent-draw null", {
  # The reason the null had to change: module medians vary several times more
  # than independent gene-set medians, so an independent null is
  # anti-conservative.
  d <- make_pool(n_bg = 600, n_block = 20, n_ctrl = 30,
                 block_lfc = 2, seed = 12)
  matched <- partner_null(d$L_ctrl, d$lfc, n_partners = 20, n_draw = 300,
                          seed = 1)
  set.seed(1)
  independent <- apply(
    matrix(sample(d$lfc, 20 * 300, replace = TRUE), nrow = 20, ncol = 300),
    2, median)
  expect_gt(sd(matched), sd(independent))
})

test_that("partner_null is reproducible under the same seed and has length n_draw", {
  d <- make_pool(n_bg = 300, n_block = 20, n_ctrl = 30,
                 block_lfc = 1, seed = 13)
  x <- partner_null(d$L_ctrl, d$lfc, n_partners = 5, n_draw = 100, seed = 7)
  y <- partner_null(d$L_ctrl, d$lfc, n_partners = 5, n_draw = 100, seed = 7)
  expect_equal(x, y)
  expect_length(x, 100)

  z <- partner_null(d$L_ctrl, d$lfc, n_partners = 5, n_draw = 100, seed = 8)
  expect_false(identical(x, z))
})

test_that("partner_null never makes a seed its own partner", {
  # If a seed could partner with itself the null would inherit the seed's own
  # log2FC, which is not part of a partner shift.
  d <- make_pool(n_bg = 300, n_block = 1, n_ctrl = 30,
                 block_lfc = 50, seed = 14)
  null <- partner_null(d$L_ctrl, d$lfc, n_partners = 20, n_draw = 300, seed = 3)
  # g1 is the only gene at lfc 50 and is uncorrelated with the rest; no draw's
  # median may be dragged anywhere near it.
  expect_lt(max(abs(null)), 1)
})

test_that("partner_p is never exactly 0 and is bounded by 1", {
  null <- rnorm(300)
  # An observed shift no null draw can reach still gets 1/(n_draw+1), not 0:
  # 300 draws cannot resolve a p below that.
  p_extreme <- partner_p(1, 1e6, null)
  expect_gt(p_extreme, 0)
  expect_equal(p_extreme, 1 / 301)
  # And the other direction is capped at 1.
  expect_equal(partner_p(1, -1e6, null), 1)
  # The negative-direction gene uses the same one-sided form.
  expect_gt(partner_p(-1, -1e6, null), 0)
})

test_that("composition_verdict boundary: share exactly 0.5 is PROGRAM, just under is MIXED", {
  expect_equal(composition_verdict(2, 1, 0.01), "PROGRAM")
  expect_equal(composition_verdict(2, 0.99999, 0.01), "MIXED")
})

test_that("composition_verdict: p >= 0.05 is always GENE-SPECIFIC regardless of share", {
  expect_equal(composition_verdict(2, 2, 0.05), "GENE-SPECIFIC")
  expect_equal(composition_verdict(2, 2, 0.2), "GENE-SPECIFIC")
})
