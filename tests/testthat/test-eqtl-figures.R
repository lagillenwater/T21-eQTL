# Data preparation behind scripts/11_eqtl_figures.R (scripts/lib/eqtl_figures.R):
# which variant each dosage panel shows, the per-gene overview rows, and the
# chr21 map bands.

library(data.table)

perm_tbl <- function() data.table(
  Gene_name = c("UP1", "DN1", "DN2"),
  n_variants = c(10L, 5L, 0L),
  best_variant = c("v_up1", "v_dn1", NA_character_),
  q_gene_bh = c(0.001, 0.4, NA_real_),
  cis_eqtl_detected = c(TRUE, FALSE, FALSE))

lane_tbl <- function() data.table(
  Gene_name = c("UP1", "DN1", "DN2", "EXP1"),
  sig_lane  = c("DE_high", "DE_low", "DE_low", "Expected_dosage"),
  eqtl_lane = c("cis_eqtl", "no_cis_eqtl", "no_GTEx_data", "not_evaluated"))

neg_tbl <- function() data.table(
  Gene_name = c("UP1", "DN1"), decoy_gene = c("DN1", "UP1"),
  n_variants = c(5L, 10L), best_variant = c("v_dn1", "v_up1"),
  q_gene_bh = c(0.9, 0.7), detected = c(FALSE, FALSE))

pos_tbl <- function() data.table(
  Gene_name = "POS1", gtex_min_p = 1e-50, n_variants = 30L,
  best_variant = "v_pos1", q_gene_bh = 0.001, detected = TRUE)

# --- panel_variants ------------------------------------------------------------

test_that("panel_variants gives each tested gene its best variant under the right group", {
  got <- panel_variants(perm_tbl(), lane_tbl(), neg_tbl(), pos_tbl())
  expect_equal(got[Gene_name == "UP1" & panel_group == "DE high", variant_id], "v_up1")
  expect_equal(got[Gene_name == "DN1" & panel_group == "DE low", variant_id], "v_dn1")
  expect_equal(got[Gene_name == "POS1", variant_id], "v_pos1")
  expect_equal(levels(got$panel_group),
               c("DE high", "DE low", "Positive control (GTEx eGenes)",
                 "Negative control (decoy variants)"))
})

test_that("panel_variants pairs a negative-control panel with the decoy set's best variant", {
  got <- panel_variants(perm_tbl(), lane_tbl(), neg_tbl(), pos_tbl())
  neg <- got[panel_group == "Negative control (decoy variants)"]
  expect_equal(neg[Gene_name == "UP1", variant_id], "v_dn1")
  expect_equal(neg[Gene_name == "UP1", decoy_gene], "DN1")
})

test_that("panel_variants drops genes with no tested variant", {
  got <- panel_variants(perm_tbl(), lane_tbl(), neg_tbl(), pos_tbl())
  expect_false("DN2" %in% got$Gene_name)
})

# --- overview_rows ---------------------------------------------------------------

test_that("overview_rows puts every deviating gene in a block with its own and decoy q", {
  got <- overview_rows(perm_tbl(), lane_tbl(), neg_tbl(), pos_tbl())
  expect_equal(as.character(got[Gene_name == "UP1", block]), "DE high")
  expect_equal(got[Gene_name == "UP1", q_own], 0.001)
  expect_equal(got[Gene_name == "UP1", q_decoy], 0.9)
  expect_equal(got[Gene_name == "DN2", eqtl_lane], "no_GTEx_data")
  expect_true(is.na(got[Gene_name == "DN2", q_own]))
  expect_false("EXP1" %in% got$Gene_name)
})

test_that("overview_rows adds positive controls as their own block without a decoy", {
  got <- overview_rows(perm_tbl(), lane_tbl(), neg_tbl(), pos_tbl())
  expect_equal(as.character(got[Gene_name == "POS1", block]), "Positive control")
  expect_equal(got[Gene_name == "POS1", q_own], 0.001)
  expect_true(is.na(got[Gene_name == "POS1", q_decoy]))
  expect_equal(levels(got$block), c("DE high", "DE low", "Positive control"))
})

# --- map_bands -------------------------------------------------------------------

test_that("map_bands places every deviating gene, taking positions from the roster first", {
  roster <- data.table(Gene_name = c("UP1", "DN1"), tss = c(10e6, 20e6))
  extra  <- data.table(Gene_name = c("DN2", "DN1"), tss = c(30e6, 99e6))
  got <- map_bands(lane_tbl(), roster, extra)
  expect_setequal(got$Gene_name, c("UP1", "DN1", "DN2"))
  expect_equal(got[Gene_name == "DN1", tss], 20e6)
  expect_equal(got[Gene_name == "DN2", tss], 30e6)
  expect_equal(got[Gene_name == "UP1", direction], "up")
  expect_equal(got[Gene_name == "DN2", eqtl_lane], "no_GTEx_data")
})

test_that("map_bands stops naming any deviating gene without a position", {
  roster <- data.table(Gene_name = "UP1", tss = 10e6)
  extra  <- data.table(Gene_name = character(0), tss = numeric(0))
  expect_error(map_bands(lane_tbl(), roster, extra), "DN1, DN2")
})

# --- best_variant_effects ----------------------------------------------------------

test_that("best_variant_effects matches lm() per panel and keeps panels separate", {
  set.seed(21)
  n <- 150
  g1 <- rbinom(n, 3, 0.3); e1 <- 0.6 * g1 + rnorm(n)
  g2 <- rbinom(n, 3, 0.4); e2 <- rnorm(n)
  long <- rbind(
    data.table(Gene_name = "A", panel_group = "DE high", variant_id = "v1", alt_dosage = g1, expr = e1),
    data.table(Gene_name = "A", panel_group = "Negative control (decoy variants)", variant_id = "v2",
               alt_dosage = g2, expr = e1))
  got <- best_variant_effects(long)
  ref <- summary(lm(e1 ~ g1))$coefficients["g1", ]
  own <- got[panel_group == "DE high"]
  expect_equal(own$slope, unname(ref["Estimate"]), tolerance = 1e-8)
  expect_equal(own$se, unname(ref["Std. Error"]), tolerance = 1e-8)
  expect_equal(own$p, unname(ref["Pr(>|t|)"]), tolerance = 1e-8)
  expect_equal(own$n, n)
  expect_equal(nrow(got), 2L)
})

test_that("best_variant_effects returns NA for a monomorphic variant", {
  long <- data.table(Gene_name = "A", panel_group = "DE low", variant_id = "v",
                     alt_dosage = rep(1, 20), expr = rnorm(20))
  got <- best_variant_effects(long)
  expect_true(is.na(got$slope))
})

# --- attach_effects ------------------------------------------------------------------

test_that("attach_effects puts the own-variant and decoy effects on each gene's row", {
  rows <- overview_rows(perm_tbl(), lane_tbl(), neg_tbl(), pos_tbl())
  eff <- data.table(
    Gene_name   = c("UP1", "UP1", "DN1", "POS1"),
    panel_group = c("DE high", "Negative control (decoy variants)", "DE low",
                    "Positive control (GTEx eGenes)"),
    variant_id  = c("v_up1", "v_dn1", "v_dn1", "v_pos1"),
    n = 100L, slope = c(-0.5, 0.1, 0.2, 1.2), se = c(0.05, 0.04, 0.06, 0.1), p = 0.01)
  got <- attach_effects(rows, eff)
  expect_equal(got[Gene_name == "UP1", own_slope], -0.5)
  expect_equal(got[Gene_name == "UP1", own_se], 0.05)
  expect_equal(got[Gene_name == "UP1", decoy_slope], 0.1)
  expect_equal(got[Gene_name == "POS1", own_slope], 1.2)
  expect_true(is.na(got[Gene_name == "POS1", decoy_slope]))
  expect_true(is.na(got[Gene_name == "DN2", own_slope]))
  expect_equal(nrow(got), nrow(rows))
})
