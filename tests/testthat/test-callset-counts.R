test_that("allele_type follows htslib's typing of REF/ALT pairs", {
  expect_equal(allele_type("A", "G"), "snp")
  expect_equal(allele_type("A", "A"), "ref")
  expect_equal(allele_type("A", "*"), "overlap")
  expect_equal(allele_type("A", "AT"), "indel")            # insertion
  expect_equal(allele_type("ATG", "A"), "indel")           # deletion
  expect_equal(allele_type("AT", "GT"), "snp")             # trims to one base each
  expect_equal(allele_type("AC", "GT"), "mnp")
  expect_equal(allele_type("ACG", "T"), "other")           # length differs, flanks do not match
  expect_equal(allele_type("CAG", "CTTG"), "other")        # substitution plus insertion
  expect_equal(allele_type("CAG", "CAAG"), "indel")
})

test_that("callset_counts counts a record under every type its alleles have", {
  got <- callset_counts(ref = c("A", "AT", "C", "G"), alt = c("G", "A", "T,CA", "GTT"))
  expect_equal(got$records, 4L)
  expect_equal(got$with_snv, 2L)     # A>G; C>T
  expect_equal(got$with_indel, 3L)   # AT>A; C>CA; G>GTT
  expect_equal(got$with_both, 1L)    # the multiallelic C record
})
