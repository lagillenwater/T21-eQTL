library(data.table)

LFC_PLOIDY <- log2(1.5)

# Eight genes: chr21 protein-coding x2 (kept), chr21 lncRNA (dropped by the
# protein-coding switch), chr22 protein-coding x3 (kept), a chr21 gene with an
# NA uncorrected estimate (dropped from BOTH scales), and a chr1 gene (not
# requested). chr21 corrected = uncorrected - log2(1.5); chr22 unchanged.
toy_res <- function() {
  data.table(
    EnsemblID   = sprintf("G%02d", 1:8),
    Chr         = c("chr21", "chr21", "chr21", "chr22", "chr22", "chr22",
                    "chr21", "chr1"),
    Gene_type   = c("protein_coding", "protein_coding", "lncRNA",
                    "protein_coding", "protein_coding", "protein_coding",
                    "protein_coding", "protein_coding"),
    raw_log2FC  = c(0.6, 0.5, 0.7, 0.1, -0.2, 0.0, NA, 0.3),
    norm_log2FC = c(0.6 - LFC_PLOIDY, 0.5 - LFC_PLOIDY, 0.7 - LFC_PLOIDY,
                    0.1, -0.2, 0.0, 0.2, 0.3)
  )
}

test_that("each retained gene appears once per scale on its own chromosome", {
  long <- ploidy_distribution_long(toy_res())
  expect_equal(nrow(long), 10)
  expect_equal(long[, .N, by = .(chromosome, scale)][order(chromosome, scale), N],
               c(2, 2, 3, 3))
  expect_equal(levels(long$chromosome), c("chr21", "chr22"))
  expect_equal(levels(long$scale), c("uncorrected", "ploidy-corrected"))
  expect_false("G07" %in% long$EnsemblID)
  expect_false("G08" %in% long$EnsemblID)
  expect_false("G03" %in% long$EnsemblID)
})

test_that("the protein-coding restriction is a switch", {
  long <- ploidy_distribution_long(toy_res(), protein_coding_only = FALSE)
  expect_true("G03" %in% long$EnsemblID)
  expect_equal(nrow(long), 12)
})

test_that("a requested chromosome with no genes is an error, not an empty curve", {
  expect_error(ploidy_distribution_long(toy_res(), chromosomes = c("chr21", "22")),
               "22")
})

test_that("input columns are checked", {
  expect_error(ploidy_distribution_long(toy_res()[, !"norm_log2FC"]),
               "norm_log2FC")
})

test_that("shift stats recover the rigid log2(1.5) shift on chr21 and none on chr22", {
  sh <- ploidy_shift_stats(ploidy_distribution_long(toy_res()))
  expect_equal(as.character(sh$chromosome), c("chr21", "chr22"))
  expect_equal(sh$n, c(2, 3))
  expect_equal(sh[chromosome == "chr21", median_shift], -LFC_PLOIDY)
  expect_equal(sh[chromosome == "chr21", c(min_shift, max_shift)],
               c(-LFC_PLOIDY, -LFC_PLOIDY))
  expect_equal(sh[chromosome == "chr22", c(median_shift, min_shift, max_shift)],
               c(0, 0, 0))
})

test_that("per-set summary stats are one row per chromosome x scale", {
  st <- ploidy_distribution_stats(ploidy_distribution_long(toy_res()))
  expect_equal(nrow(st), 4)
  expect_equal(st$n, c(2, 2, 3, 3))
  expect_equal(st[chromosome == "chr21" & scale == "uncorrected", median], 0.55)
  expect_equal(st[chromosome == "chr21" & scale == "ploidy-corrected", median],
               0.55 - LFC_PLOIDY)
  expect_equal(st[chromosome == "chr22" & scale == "uncorrected", median],
               st[chromosome == "chr22" & scale == "ploidy-corrected", median])
})
