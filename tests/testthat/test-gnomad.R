# Tests for scripts/lib/gnomad.R. Only the pure parts are exercised: the
# region planner, the INFO parser and the cache reader. Nothing here touches
# the network - fetch_gnomad_af() is the single function that does, and its
# inputs and outputs are these three.

library(data.table)

vcf_line <- function(pos, ref, alt, info, filter = "PASS")
  paste("chr21", pos, ".", ref, alt, ".", filter, info, sep = "\t")

# --- gnomad_regions ----------------------------------------------------------

test_that("gnomad_regions merges positions within the gap and splits beyond it", {
  got <- gnomad_regions(c(1000, 1200, 6000, 6001, 50000), gap = 2000)
  expect_equal(nrow(got), 3L)
  expect_equal(got$start, c(1000, 6000, 50000))
  expect_equal(got$end,   c(1200, 6001, 50000))
})

test_that("gnomad_regions is one range per position at gap 0", {
  expect_equal(nrow(gnomad_regions(c(10, 12, 14), gap = 0)), 3L)
})

test_that("gnomad_regions drops NA, de-duplicates and handles no positions", {
  got <- gnomad_regions(c(500, NA, 500), gap = 100)
  expect_equal(nrow(got), 1L)
  expect_equal(nrow(gnomad_regions(integer(0))), 0L)
  expect_equal(nrow(gnomad_regions(NA_integer_)), 0L)
})

# --- gnomad_info_field -------------------------------------------------------

test_that("gnomad_info_field reads a key without matching a longer key", {
  info <- "AC=9;AN=100;AF=0.09;AF_nfe=0.12;AF_XX=0.4;nhomalt=1"
  expect_equal(gnomad_info_field(info, "AF"), "0.09")
  expect_equal(gnomad_info_field(info, "AF_nfe"), "0.12")
  expect_equal(gnomad_info_field("AF=0.5;AC=2", "AF"), "0.5")   # first field
})

test_that("gnomad_info_field returns NA for a key that is absent", {
  expect_true(is.na(gnomad_info_field("AC=9;AN=100", "AF")))
})

# --- parse_gnomad_lines ------------------------------------------------------

test_that("parse_gnomad_lines returns one row per line with numeric frequencies", {
  got <- parse_gnomad_lines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    vcf_line(13858796, "C", "G", "AC=9665;AN=152094;AF=0.0635462;AF_nfe=0.0212"),
    vcf_line(13858800, "A", "T", "AC=1;AN=100;AF=0.6;AF_nfe=0.55", filter = "AC0")))
  expect_equal(nrow(got), 2L)
  expect_equal(got$POS, c(13858796L, 13858800L))
  expect_equal(got$REF, c("C", "A"))
  expect_equal(got$AF, c(0.0635462, 0.6))
  expect_equal(got$AF_nfe, c(0.0212, 0.55))
  expect_equal(got$filter, c("PASS", "AC0"))
})

test_that("parse_gnomad_lines drops a multi-allelic line rather than guessing", {
  got <- parse_gnomad_lines(vcf_line(100, "A", "T,G", "AF=0.2"))
  expect_equal(nrow(got), 0L)
})

test_that("parse_gnomad_lines gives NA where the frequency key is missing", {
  got <- parse_gnomad_lines(vcf_line(100, "A", "T", "AC=3;AN=90"))
  expect_true(is.na(got$AF))
})

test_that("parse_gnomad_lines returns a typed empty table for no lines", {
  got <- parse_gnomad_lines(character(0))
  expect_equal(nrow(got), 0L)
  expect_true(all(c("POS", "REF", "ALT", "filter", GNOMAD_AF_FIELDS) %in% names(got)))
})

# --- read_gnomad_cache -------------------------------------------------------

test_that("read_gnomad_cache reads the rows and the queried positions", {
  dir <- withr::local_tempdir()
  af  <- file.path(dir, "af.csv"); pos <- file.path(dir, "pos.txt")
  fwrite(data.table(POS = 100L, REF = "A", ALT = "T", filter = "PASS",
                    AF = 0.1, AF_nfe = 0.2), af)
  writeLines(c("100", "200"), pos)
  got <- read_gnomad_cache(af, pos)
  expect_equal(nrow(got$af), 1L)
  expect_equal(got$queried, c(100L, 200L))
})

test_that("read_gnomad_cache on an absent cache is empty, not an error", {
  got <- read_gnomad_cache(tempfile(), tempfile())
  expect_equal(nrow(got$af), 0L)
  expect_equal(length(got$queried), 0L)
})

test_that("ensure_gnomad_af stays offline and reports the positions it lacks", {
  dir <- withr::local_tempdir()
  af  <- file.path(dir, "af.csv"); pos <- file.path(dir, "pos.txt")
  fwrite(data.table(POS = 100L, REF = "A", ALT = "T", filter = "PASS",
                    AF = 0.1, AF_nfe = 0.2), af)
  writeLines("100", pos)
  withr::local_envvar(T21_GNOMAD_OFFLINE = "1")
  got <- ensure_gnomad_af(c(100L, 200L), af, pos, verbose = FALSE)
  expect_equal(got$POS, 100L)
  expect_match(attr(got, "gnomad_status"), "1 of 2 positions not queried")
})

test_that("ensure_gnomad_af reports ok when every position was queried", {
  dir <- withr::local_tempdir()
  af  <- file.path(dir, "af.csv"); pos <- file.path(dir, "pos.txt")
  fwrite(data.table(POS = 100L, REF = "A", ALT = "T", filter = "PASS",
                    AF = 0.1, AF_nfe = 0.2), af)
  writeLines(c("100", "200"), pos)   # 200 queried, gnomAD carries no variant there
  withr::local_envvar(T21_GNOMAD_OFFLINE = "1")
  got <- ensure_gnomad_af(c(100L, 200L), af, pos, verbose = FALSE)
  expect_equal(attr(got, "gnomad_status"), "ok")
})
