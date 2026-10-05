# callset_counts.R
#
# Variant-type counts of a genotype call set, made the way `bcftools stats`
# makes them, so the Methods counts can be checked against the files the
# analysis reads without bcftools. Each ALT allele is typed by htslib's
# bcf_set_variant_type rule; a record is counted under every type any of its
# alleles has, so the SNV and indel counts of a multiallelic call set overlap.

#' Type one REF/ALT allele pair as htslib does: strip the shared prefix; if
#' one allele is used up the pair is an indel; otherwise strip the shared
#' suffix, and one base against one base is an SNP, one base against several
#' an indel when the flanking bases match (else other), equal lengths an MNP,
#' anything else other. "*" (spanning deletion) is an overlap.
#' @return one of "snp", "mnp", "indel", "other", "overlap", "ref"
allele_type <- function(ref, alt) {
  if (alt == "*") return("overlap")
  if (nchar(ref) == 1 && nchar(alt) == 1) return(if (toupper(ref) == toupper(alt)) "ref" else "snp")
  r <- strsplit(toupper(ref), "")[[1]]
  a <- strsplit(toupper(alt), "")[[1]]
  i <- 1L
  while (i <= length(r) && i <= length(a) && r[i] == a[i]) i <- i + 1L
  if (i > length(r) && i > length(a)) return("ref")
  if (i > length(r) || i > length(a)) return("indel")
  re <- length(r); ae <- length(a)
  while (re > i && ae > i && r[re] == a[ae]) { re <- re - 1L; ae <- ae - 1L }
  if (ae == i) {
    if (re == i) return("snp")
    return(if (r[re] == a[ae]) "indel" else "other")
  }
  if (re == i) return(if (r[re] == a[ae]) "indel" else "other")
  if (re - i == ae - i) "mnp" else "other"
}

#' Record counts by type for a table of REF and ALT (comma-separated ALTs).
#' Each distinct allele pair is typed once.
#' @return data.table(records, with_snv, with_indel, with_both)
callset_counts <- function(ref, alt) {
  long <- data.table::data.table(rec = rep(seq_along(alt), lengths(strsplit(alt, ",", fixed = TRUE))),
                                 REF = rep(ref, lengths(strsplit(alt, ",", fixed = TRUE))),
                                 ALT = unlist(strsplit(alt, ",", fixed = TRUE)))
  pairs <- unique(long[, .(REF, ALT)])
  pairs[, type := mapply(allele_type, REF, ALT, USE.NAMES = FALSE)]
  long <- merge(long, pairs, by = c("REF", "ALT"))
  per <- long[, .(snv = any(type == "snp"), indel = any(type == "indel")), by = rec]
  data.table::data.table(records = length(alt), with_snv = sum(per$snv),
                         with_indel = sum(per$indel), with_both = sum(per$snv & per$indel))
}
