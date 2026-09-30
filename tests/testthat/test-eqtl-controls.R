# Standalone controls for the within-T21 cis-eQTL test (scripts/lib/eqtl_controls.R):
# TSS lookup from GTEx allpairs rows, positive-control gene selection, decoy
# (unlinked) variant-set pairing, and the shared gene-level permutation runner
# that the observed, negative and positive sets all go through.

library(data.table)

# --- gene_tss -----------------------------------------------------------------

test_that("gene_tss recovers the TSS as position minus tss_distance, per gene", {
  ap <- data.table(ensembl_stable = c("g1", "g1", "g2"),
                   POS = c(1000L, 1500L, 9000L),
                   tss_distance = c(-200L, 300L, 1000L))
  got <- gene_tss(ap)
  expect_equal(got[ensembl_stable == "g1", tss], 1200)
  expect_equal(got[ensembl_stable == "g2", tss], 8000)
})

# --- select_positive_controls -----------------------------------------------

test_that("select_positive_controls ranks eligible genes by GTEx min p and drops excluded ids", {
  candidates <- data.table(ensembl_stable = c("a", "b", "c", "d"),
                           Gene_name = c("A", "B", "C", "D"))
  gtex_min_p <- data.table(ensembl_stable = c("a", "b", "c", "d"),
                           gtex_min_p = c(1e-3, 1e-30, 1e-10, 1e-20))
  got <- select_positive_controls(candidates, gtex_min_p, exclude = "b", n = 2)
  expect_equal(got$Gene_name, c("D", "C"))
  expect_equal(got$gtex_min_p, c(1e-20, 1e-10))
})

test_that("select_positive_controls ignores candidates without GTEx coverage and caps at n", {
  candidates <- data.table(ensembl_stable = c("a", "b", "c"),
                           Gene_name = c("A", "B", "C"))
  gtex_min_p <- data.table(ensembl_stable = c("a", "c"), gtex_min_p = c(1e-5, 1e-8))
  got <- select_positive_controls(candidates, gtex_min_p, exclude = character(0), n = 5)
  expect_equal(got$Gene_name, c("C", "A"))
})

# --- assign_decoys ------------------------------------------------------------

test_that("assign_decoys picks, among genes at least min_distance away, the one with the closest variant count", {
  tss  <- c(G1 = 10e6, G2 = 12e6, G3 = 30e6, G4 = 40e6)
  nvar <- c(G1 = 100,  G2 = 90,   G3 = 120,  G4 = 1)
  got <- assign_decoys(tss, nvar, min_distance = 5e6)
  # G2 is too close to G1; of G3 (120) and G4 (1), G3 is the closer size
  expect_equal(got[Gene_name == "G1", decoy_gene], "G3")
  expect_equal(got[Gene_name == "G1", distance], 20e6)
  # G4's candidates are all far enough; G2 (90) is the closest size to 1
  expect_equal(got[Gene_name == "G4", decoy_gene], "G2")
})

test_that("assign_decoys breaks a variant-count tie by distance, farthest first", {
  tss  <- c(G1 = 10e6, G3 = 30e6, G4 = 40e6)
  nvar <- c(G1 = 50, G3 = 50, G4 = 50)
  got <- assign_decoys(tss, nvar, min_distance = 5e6)
  expect_equal(got[Gene_name == "G1", decoy_gene], "G4")
})

test_that("assign_decoys gives NA when no other gene is far enough", {
  tss  <- c(G1 = 10e6, G2 = 11e6)
  got <- assign_decoys(tss, c(G1 = 5, G2 = 5), min_distance = 5e6)
  expect_true(all(is.na(got$decoy_gene)))
  expect_equal(nrow(got), 2L)
})

# --- gene_level_tests -----------------------------------------------------------

sim_genotypes <- function(n, m, seed) {
  set.seed(seed)
  G <- matrix(rbinom(n * m, 3, 0.3), nrow = n)
  colnames(G) <- paste0("v", seq_len(m))
  G
}

test_that("gene_level_tests detects a spiked cis effect and not a null gene", {
  G <- sim_genotypes(n = 200, m = 6, seed = 1)
  set.seed(2)
  e_signal <- 0.8 * G[, "v3"] + rnorm(200)
  e_null   <- rnorm(200)
  expr_of  <- function(g) if (g == "SIG") e_signal else e_null
  variants_of <- list(SIG = colnames(G), NUL = colnames(G))
  got <- gene_level_tests(c("SIG", "NUL"), variants_of, G, expr_of,
                          n_perm = 200, seed_base = 100L, fdr = 0.05)
  expect_equal(got$Gene_name, c("SIG", "NUL"))
  expect_equal(got[Gene_name == "SIG", best_variant], "v3")
  expect_true(got[Gene_name == "SIG", detected])
  expect_false(got[Gene_name == "NUL", detected])
  expect_equal(got$n_variants, c(6L, 6L))
})

test_that("gene_level_tests emits an all-NA row for a gene with no usable variants", {
  G <- sim_genotypes(n = 50, m = 2, seed = 3)
  got <- gene_level_tests("X", list(X = character(0)), G, function(g) rnorm(50),
                          n_perm = 50, seed_base = 1L, fdr = 0.05)
  expect_equal(got$n_variants, 0L)
  expect_true(is.na(got$p_gene_perm))
  expect_false(got$detected)
})

test_that("gene_level_tests is reproducible for a fixed seed_base", {
  G <- sim_genotypes(n = 80, m = 4, seed = 4)
  set.seed(5); e <- rnorm(80)
  a <- gene_level_tests("A", list(A = colnames(G)), G, function(g) e,
                        n_perm = 100, seed_base = 7L, fdr = 0.05)
  b <- gene_level_tests("A", list(A = colnames(G)), G, function(g) e,
                        n_perm = 100, seed_base = 7L, fdr = 0.05)
  expect_equal(a$p_gene_perm, b$p_gene_perm)
})

# --- match_positive_controls (controls v2) ------------------------------------

test_that("match_positive_controls picks the nearest candidate in effect-size and expression", {
  targets <- data.table(Gene_name = c("T1", "T2"), abs_afc = c(1.0, 0.3), baseMean = c(100, 2000))
  candidates <- data.table(ensembl_stable = c("a", "b", "c", "d"), Gene_name = c("A", "B", "C", "D"),
                           abs_afc = c(1.1, 0.3, 4, 0.9), baseMean = c(120, 2100, 100, 5000),
                           tss = c(1e6, 10e6, 20e6, 30e6))
  got <- match_positive_controls(targets, candidates, min_separation = 1e6)
  expect_equal(got$target_gene, c("T1", "T2"))
  expect_equal(got$Gene_name, c("A", "B"))
})

test_that("match_positive_controls uses one locus once and skips NA targets", {
  targets <- data.table(Gene_name = c("T1", "T2", "T3"), abs_afc = c(1, 1, NA), baseMean = c(100, 100, 100))
  candidates <- data.table(ensembl_stable = c("a", "b", "c"), Gene_name = c("A", "B", "C"),
                           abs_afc = c(1, 1, 0.5), baseMean = c(100, 100, 100),
                           tss = c(5e6, 5.4e6, 20e6))   # A and B share a locus
  got <- match_positive_controls(targets, candidates, min_separation = 1e6)
  expect_equal(nrow(got), 2)
  expect_setequal(got$Gene_name, c("A", "C"))          # B is within 1 Mb of A
})

# --- assign_decoy_sets --------------------------------------------------------

test_that("assign_decoy_sets returns k distant donors per gene, closest variant count first", {
  tss <- c(G1 = 10e6, G2 = 40e6); nv <- c(G1 = 100, G2 = 50)
  dtss <- c(G1 = 10e6, G2 = 40e6, P1 = 12e6, P2 = 30e6, P3 = 45e6, P4 = 20e6)
  dnv  <- c(G1 = 100, G2 = 50, P1 = 100, P2 = 90, P3 = 55, P4 = 100)
  got <- assign_decoy_sets(tss, nv, dtss, dnv, min_distance = 5e6, k = 2)
  expect_equal(got[Gene_name == "G1", decoy_gene], c("P4", "P2"))   # P1 too close, G1 itself excluded
  expect_equal(got[Gene_name == "G2", decoy_gene], c("P3", "P2"))
  expect_equal(got[Gene_name == "G1", decoy_rank], 1:2)
  expect_true(all(got$distance >= 5e6))
})

test_that("assign_decoy_sets returns fewer than k when few donors are far enough", {
  tss <- c(G1 = 10e6); nv <- c(G1 = 10)
  got <- assign_decoy_sets(tss, nv, c(G1 = 10e6, P1 = 12e6, P2 = 30e6), c(G1 = 10, P1 = 10, P2 = 10),
                           min_distance = 5e6, k = 3)
  expect_equal(got$decoy_gene, "P2")
})
