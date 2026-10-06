# Data preparation behind scripts/11_eqtl_figures.R (scripts/lib/eqtl_figures.R):
# which variant each dosage panel shows, which allele it is drawn on, and the
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

# --- panel_variants ------------------------------------------------------------

test_that("panel_variants gives each tested deviating gene its best variant under its block", {
  got <- panel_variants(perm_tbl(), lane_tbl())
  expect_equal(got[Gene_name == "UP1", variant_id], "v_up1")
  expect_equal(as.character(got[Gene_name == "UP1", panel_group]), "DE high")
  expect_equal(got[Gene_name == "DN1", variant_id], "v_dn1")
  expect_equal(as.character(got[Gene_name == "DN1", panel_group]), "DE low")
  expect_equal(levels(got$panel_group), c("DE high", "DE low"))
})

test_that("panel_variants drops genes with no tested variant", {
  got <- panel_variants(perm_tbl(), lane_tbl())
  expect_false("DN2" %in% got$Gene_name)
  expect_equal(nrow(got), 2L)
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
    data.table(Gene_name = "A", panel_group = "DE high", variant_id = "v2", alt_dosage = g2, expr = e1))
  got <- best_variant_effects(long)
  ref <- summary(lm(e1 ~ g1))$coefficients["g1", ]
  own <- got[variant_id == "v1"]
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

# --- orient_panels -------------------------------------------------------------
# A deviating-gene panel is drawn on the allele whose GTEx effect runs the same
# way as the gene's deviation, so the panel trends with its block label. That
# allele is the minor one exactly when the minor allele runs with the
# deviation, which is what makes plot_allele_role double as the badge.

variant_tbl <- function() data.table(
  variant_id   = c("v_up1", "v_dn1"),
  minor_allele = c("T", "A"),
  major_allele = c("C", "G"),
  alt_is_minor = c(TRUE, FALSE),
  gtex_maf     = c(0.29, 0.13))

gtex_tbl <- function() data.table(
  variant_id       = c("v_up1", "v_dn1"),
  Gene_name        = c("UP1", "DN1"),
  gtex_slope_minor = c(0.3, 0.25))

gene_dir_tbl <- function() data.table(
  Gene_name = c("UP1", "DN1", "DN2"), deviation_sign = c(1, -1, -1))

oriented <- function() orient_panels(panel_variants(perm_tbl(), lane_tbl()),
                                     variant_tbl(), gtex_tbl(), gene_dir_tbl())

test_that("a DE high gene whose minor allele raises expression is drawn on the minor allele", {
  got <- oriented()[Gene_name == "UP1" & panel_group == "DE high"]
  expect_true(got$with_deviation)
  expect_true(got$plot_minor)
  expect_equal(got$plot_allele, "T")
  expect_equal(got$plot_allele_role, "minor")
  expect_equal(got$canonical_dir, "raises")
})

test_that("a DE low gene whose minor allele raises expression is drawn on the major allele", {
  got <- oriented()[Gene_name == "DN1" & panel_group == "DE low"]
  expect_false(got$with_deviation)     # minor A raises, the gene deviates low
  expect_false(got$plot_minor)
  expect_equal(got$plot_allele, "G")
  expect_equal(got$plot_allele_role, "major")
  expect_equal(got$canonical_dir, "raises")
})

test_that("a panel with no GTEx slope for its gene has nothing to match and stays on the minor allele", {
  got <- orient_panels(panel_variants(perm_tbl(), lane_tbl()), variant_tbl(),
                       gtex_tbl()[Gene_name != "DN1"], gene_dir_tbl())[Gene_name == "DN1"]
  expect_true(is.na(got$canonical_dir))
  expect_true(is.na(got$with_deviation))
  expect_true(got$plot_minor)
  expect_equal(got$plot_allele, "A")
})

test_that("orientation follows GTEx, not the sign of the data being plotted", {
  # Same panel, GTEx slope flipped: the plotted allele has to flip with it.
  g <- gtex_tbl(); g[Gene_name == "UP1", gtex_slope_minor := -0.3]
  got <- orient_panels(panel_variants(perm_tbl(), lane_tbl()),
                       variant_tbl(), g, gene_dir_tbl())[Gene_name == "UP1" & panel_group == "DE high"]
  expect_false(got$plot_minor)
  expect_equal(got$plot_allele, "C")
})

test_that("orient_panels keeps one row per panel", {
  p <- panel_variants(perm_tbl(), lane_tbl())
  expect_equal(nrow(oriented()), nrow(p))
})

# --- panel_dosage --------------------------------------------------------------

test_that("panel_dosage is minor dosage when the panel shows the minor allele", {
  expect_equal(panel_dosage(0:3, TRUE, TRUE, 3), as.numeric(0:3))
  expect_equal(panel_dosage(0:3, FALSE, TRUE, 3), c(3, 2, 1, 0))
})

test_that("panel_dosage is the major-allele count when the panel shows the major allele", {
  expect_equal(panel_dosage(0:3, TRUE, FALSE, 3), c(3, 2, 1, 0))
  expect_equal(panel_dosage(0:3, FALSE, FALSE, 3), as.numeric(0:3))
})

test_that("panel_dosage recycles a scalar orientation over a dosage vector", {
  expect_equal(panel_dosage(c(0, 1, 2, 3), TRUE, FALSE, 3), c(3, 2, 1, 0))
  expect_true(all(is.na(panel_dosage(0:3, NA, TRUE, 3))))
})

# --- stack_bands -----------------------------------------------------------------

test_that("stack_bands gives a lone band the full height", {
  got <- stack_bands(c(10, 20, 30), half_width = 0.1)
  expect_equal(got$ymin, c(0, 0, 0))
  expect_equal(got$ymax, c(1, 1, 1))
})

test_that("stack_bands splits overlapping bands into equal slices, in position order", {
  got <- stack_bands(c(33.07, 33.03, 33.071, 40), half_width = 0.12)
  # 33.03, 33.07, 33.071 overlap: bottom to top in position order
  expect_equal(got$ymin, c(1, 0, 2, 0) / c(3, 3, 3, 1))
  expect_equal(got$ymax, c(2, 1, 3, 1) / c(3, 3, 3, 1))
})

test_that("stack_bands chains overlaps into one cluster", {
  got <- stack_bands(c(1.0, 1.2, 1.4), half_width = 0.11)   # 1.0~1.2 and 1.2~1.4 overlap; 1.0 and 1.4 do not
  expect_equal(got$ymax - got$ymin, rep(1 / 3, 3))
})

# --- map_outcome -----------------------------------------------------------------

test_that("map_outcome splits the genes without a variant to test by GTEx status", {
  got <- map_outcome(c("cis_eqtl", "no_cis_eqtl", "no_GTEx_data", "no_GTEx_data"),
                     gtex_tested = c(TRUE, TRUE, TRUE, FALSE))
  expect_equal(as.character(got), c("cis-eQTL detected", "tested, not detected",
                                    "GTEx tested, no eQTL", "not tested by GTEx"))
  expect_equal(levels(got), MAP_OUTCOMES)
})

test_that("map_outcome refuses an unknown eQTL lane", {
  expect_error(map_outcome("maybe", TRUE), "maybe")
})
