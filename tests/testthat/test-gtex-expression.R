# Tests for scripts/lib/gtex_expression.R. Only the pure parts are exercised:
# the GCT reader and the per-gene tissue summary. Nothing here touches the
# network - gtex_median_tpm_file() is the single function that does.

library(data.table)

write_gct <- function(rows, tissues = c("Liver", "Whole_Blood", "Testis"), gz = FALSE) {
  path <- tempfile(fileext = if (gz) ".gct.gz" else ".gct")
  con <- if (gz) gzfile(path, "w") else file(path, "w")
  writeLines(c("#1.2", sprintf("%d\t%d", length(rows), length(tissues)),
               paste(c("Name", "Description", tissues), collapse = "\t"), rows), con)
  close(con)
  path
}

# --- read_gtex_median_tpm ----------------------------------------------------

test_that("read_gtex_median_tpm strips the Ensembl version and keeps the tissue columns", {
  got <- read_gtex_median_tpm(write_gct(c("ENSG00000159216.19\tRUNX1\t1.5\t6.9\t0.2",
                                          "ENSG00000142173.15\tCOL6A2\t10\t8.1\t3")))
  expect_equal(got$ensembl_stable, c("ENSG00000159216", "ENSG00000142173"))
  expect_equal(got$gene_name, c("RUNX1", "COL6A2"))
  expect_equal(got$Whole_Blood, c(6.9, 8.1))
  expect_setequal(setdiff(names(got), c("ensembl_stable", "gene_name")), c("Liver", "Whole_Blood", "Testis"))
})

test_that("read_gtex_median_tpm drops the PAR_Y copies so each stable id appears once", {
  got <- read_gtex_median_tpm(write_gct(c("ENSG00000002586.20\tCD99\t50\t60\t70",
                                          "ENSG00000002586.20_PAR_Y\tCD99\t0\t0\t0")))
  expect_equal(got$ensembl_stable, "ENSG00000002586")
  expect_equal(got$Whole_Blood, 60)
})

test_that("read_gtex_median_tpm reads a gzipped file", {
  got <- read_gtex_median_tpm(write_gct("ENSG00000159216.19\tRUNX1\t1.5\t6.9\t0.2", gz = TRUE))
  expect_equal(got$Whole_Blood, 6.9)
})

test_that("read_gtex_median_tpm refuses a file without a Whole_Blood column", {
  path <- write_gct("ENSG00000159216.19\tRUNX1\t1.5\t0.2", tissues = c("Liver", "Testis"))
  expect_error(read_gtex_median_tpm(path), "Whole_Blood")
})

test_that("read_gtex_median_tpm refuses a file that is not GCT", {
  path <- tempfile(fileext = ".gct")
  writeLines(c("<html>Access denied</html>"), path)
  expect_error(read_gtex_median_tpm(path), "GCT")
})

# --- tissue_expression_summary -----------------------------------------------

tpm <- data.table(ensembl_stable = c("G1", "G2", "G3", "G4"),
                  Liver       = c(10, 1, 0, 4),
                  Whole_Blood = c( 1, 1, 0, 4),
                  Testis      = c( 0, 1, 0, 4),
                  Ovary       = c( 5, 1, 0, 8))

test_that("tissue_expression_summary gives the focal TPM and its rank among tissues", {
  got <- tissue_expression_summary(tpm, focal = "Whole_Blood")
  expect_equal(got$ensembl_stable, tpm$ensembl_stable)
  expect_equal(got$focal_tpm, c(1, 1, 0, 4))
  expect_equal(got$n_tissues, rep(4L, 4))
  expect_equal(got$focal_rank[1], 3)          # 10 and 5 above it
})

test_that("tissue_expression_summary averages the rank over ties", {
  got <- tissue_expression_summary(tpm, focal = "Whole_Blood")
  expect_equal(got$focal_rank[2], 2.5)        # four-way tie at 1
  expect_equal(got$focal_rank[4], 3)          # 8 above, three-way tie at 4: ranks 2-4
})

test_that("tissue_expression_summary leaves rank, top tissue and tau NA for a gene expressed nowhere", {
  got <- tissue_expression_summary(tpm, focal = "Whole_Blood")
  expect_true(is.na(got$focal_rank[3]))
  expect_true(is.na(got$top_tissue[3]))
  expect_true(is.na(got$tau[3]))
  expect_equal(got$top_tpm[3], 0)
})

test_that("tissue_expression_summary names the top tissue and the focal fraction of it", {
  got <- tissue_expression_summary(tpm, focal = "Whole_Blood")
  expect_equal(got$top_tissue[c(1, 4)], c("Liver", "Ovary"))
  expect_equal(got$top_tpm[c(1, 4)], c(10, 8))
  expect_equal(got$focal_frac_of_top[c(1, 4)], c(0.1, 0.5))
})

test_that("tissue_expression_summary computes tau on log2(TPM + 1)", {
  got <- tissue_expression_summary(tpm, focal = "Whole_Blood")
  x <- log2(c(10, 1, 0, 5) + 1)
  expect_equal(got$tau[1], sum(1 - x / max(x)) / (length(x) - 1))
  expect_equal(got$tau[2], 0)                 # uniform: not specific
  one_tissue <- data.table(ensembl_stable = "G", A = 0, B = 0, Whole_Blood = 7)
  expect_equal(tissue_expression_summary(one_tissue)$tau, 1)
})

test_that("tissue_expression_summary stops when the focal tissue is missing", {
  expect_error(tissue_expression_summary(tpm, focal = "Spleen"), "Spleen")
})

# --- read_gtf_exons ------------------------------------------------------------

gtf_line <- function(feature, start, end, gene_id)
  paste("chr21", "HAVANA", feature, start, end, ".", "+", ".",
        sprintf('gene_id "%s"; transcript_id "ENST1.1"; gene_name "X";', gene_id), sep = "\t")

test_that("read_gtf_exons keeps exon rows only, with the stable gene id", {
  path <- tempfile(fileext = ".gtf.gz")
  con <- gzfile(path, "w")
  writeLines(c("##description: test", gtf_line("gene", 100, 900, "ENSG00000159216.19"),
               gtf_line("exon", 100, 200, "ENSG00000159216.19"),
               gtf_line("exon", 500, 900, "ENSG00000159216.19"),
               gtf_line("exon", 10, 20, "ENSG00000002586.20_PAR_Y")), con)
  close(con)
  got <- read_gtf_exons(path)
  expect_equal(got$ensembl_stable, rep("ENSG00000159216", 2))
  expect_equal(got$start, c(100L, 500L))
  expect_equal(got$end, c(200L, 900L))
})

# --- exon_union_length ---------------------------------------------------------

test_that("exon_union_length counts overlapping exons once and sums disjoint ones", {
  ex <- data.table(ensembl_stable = c("A", "A", "A", "B", "B"),
                   start = c(100L, 150L, 500L, 10L, 30L),
                   end   = c(200L, 250L, 599L, 40L, 35L))   # B: second exon nested in the first
  got <- exon_union_length(ex)
  expect_equal(got[ensembl_stable == "A", length_bp], 151 + 100)   # 100-250 and 500-599
  expect_equal(got[ensembl_stable == "B", length_bp], 31)
})

test_that("exon_union_length is exact for adjacent exons (no off-by-one)", {
  ex <- data.table(ensembl_stable = c("A", "A"), start = c(1L, 11L), end = c(10L, 20L))
  expect_equal(exon_union_length(ex)$length_bp, 20)
})

test_that("exon_union_length does not depend on row order", {
  ex <- data.table(ensembl_stable = c("A", "B", "A"), start = c(500L, 10L, 100L), end = c(599L, 40L, 200L))
  got <- exon_union_length(ex)
  expect_equal(got[ensembl_stable == "A", length_bp], 201)
  expect_equal(got[ensembl_stable == "B", length_bp], 31)
})

# --- counts_to_tpm -------------------------------------------------------------

test_that("counts_to_tpm gives each sample a total of one million", {
  counts <- matrix(c(10, 20, 30, 5, 0, 95), nrow = 3, dimnames = list(c("g1", "g2", "g3"), c("s1", "s2")))
  got <- counts_to_tpm(counts, c(g1 = 1000, g2 = 2000, g3 = 500))
  expect_equal(unname(colSums(got)), c(1e6, 1e6))
})

test_that("counts_to_tpm scales by length: equal counts, double length, half the TPM", {
  counts <- matrix(c(100, 100), nrow = 2, dimnames = list(c("g1", "g2"), "s1"))
  got <- counts_to_tpm(counts, c(g1 = 1000, g2 = 2000))
  expect_equal(got["g1", "s1"] / got["g2", "s1"], 2)
})

test_that("counts_to_tpm drops genes without a length", {
  counts <- matrix(c(10, 20, 30), nrow = 3, dimnames = list(c("g1", "g2", "g3"), "s1"))
  got <- counts_to_tpm(counts, c(g1 = 1000, g3 = 1000))
  expect_equal(rownames(got), c("g1", "g3"))
  expect_equal(got["g1", "s1"], 0.25e6)
})

# --- gtex_tissue_group ---------------------------------------------------------

test_that("gtex_tissue_group puts whole blood in a group of its own", {
  expect_equal(as.character(gtex_tissue_group("Whole_Blood")), "Whole blood")
})

test_that("gtex_tissue_group assigns tissues by organ system", {
  got <- as.character(gtex_tissue_group(c(
    "Brain_Cortex", "Nerve_Tibial", "Pituitary", "Spleen", "Cells_EBV-transformed_lymphocytes",
    "Heart_Left_Ventricle", "Artery_Aorta", "Muscle_Skeletal", "Colon_Transverse_Mucosa", "Liver_Hepatocyte",
    "Minor_Salivary_Gland", "Testis", "Cervix_Ectocervix", "Lung", "Cells_Cultured_fibroblasts",
    "Breast_Mammary_Tissue")))
  expect_equal(got, c(rep("Brain and nerve", 3), rep("Immune and lymphoid", 2),
                      rep("Heart, vessels and muscle", 3), rep("Digestive", 3),
                      rep("Reproductive", 2), rep("Other", 3)))
})

test_that("gtex_tissue_group sends the ileal lymphoid aggregate to immune, not digestive", {
  got <- as.character(gtex_tissue_group(c("Small_Intestine_Terminal_Ileum_Lymphode_Aggregate",
                                          "Small_Intestine_Terminal_Ileum")))
  expect_equal(got, c("Immune and lymphoid", "Digestive"))
})

test_that("gtex_tissue_group returns a factor in the fixed group order", {
  got <- gtex_tissue_group(c("Lung", "Whole_Blood"))
  expect_s3_class(got, "factor")
  expect_equal(levels(got), GTEX_TISSUE_GROUPS)
})
