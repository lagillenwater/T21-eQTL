# alleles.R
#
# Minor-allele referencing for the cis-eQTL stage.
#
# The pipeline codes genotype as ALT-allele dosage (script 02) and GTEx's
# `slope` is ALT-referenced too, so the two are internally consistent - but
# the ALT allele of a variant id is whichever base differs from the reference
# assembly, not the rarer allele in the population. Roughly a quarter of the
# retained cis variants have ALT as the MAJOR allele, so "slope per alt
# allele" points one way for some variants and the other way for others, and
# a reported direction cannot be compared with a published eQTL direction.
#
# These helpers re-express dosage and slopes per copy of the MINOR allele in a
# named reference population, which is the allele eQTL results are normally
# reported against. Nothing here changes a test statistic: flipping the sign
# of a regressor flips the sign of its slope and leaves |t|, p and R-squared
# untouched, so gene-level permutation p, q and lane calls are invariant.
#
# Allele frequencies come from the ALT-allele frequency of a reference
# population: GTEx whole blood (`af`, the allpairs column, the population the
# cis-eQTL itself was called in) as primary, gnomAD v4.1 genomes
# (scripts/lib/gnomad.R) as an independent check on which allele is minor.

# A variant whose ALT frequency is exactly 0.5 has no minor allele. The tie is
# broken toward ALT so the coding stays deterministic; `af_is_tie()` flags the
# case for reporting.
MAF_TIE_TOL <- 1e-12

# chr21 copy number in a trisomic subject. The HTP T21 calls are ALT-allele
# counts on a 0-3 scale, so this is both the denominator of an allele
# frequency estimated from those dosages and the reflection point that turns
# ALT dosage into minor-allele dosage.
T21_CHR21_PLOIDY <- 3L

#' Is the ALT allele the minor allele at this ALT allele frequency?
#' @param af numeric ALT allele frequency in [0, 1]; NA propagates.
#' @return logical, TRUE when ALT is the rarer allele (ties go to ALT).
alt_is_minor <- function(af) {
  bad <- !is.na(af) & (af < 0 | af > 1)
  if (any(bad)) stop("allele frequency outside [0, 1]: ", af[which(bad)[1]], call. = FALSE)
  ifelse(is.na(af), NA, af <= 0.5)
}

#' ALT frequency within MAF_TIE_TOL of 0.5, where "minor" is undefined.
af_is_tie <- function(af) !is.na(af) & abs(af - 0.5) <= MAF_TIE_TOL

#' Minor allele frequency from an ALT allele frequency.
maf_from_af <- function(af) pmin(af, 1 - af)

#' Recycle a set of vectors to a common length.
#'
#' base::ifelse sizes its result from the test alone, so a nested call with a
#' scalar branch silently collapses a vector to one value. Every function
#' below recycles first and only then branches.
recycle_common <- function(...) {
  v <- list(...)
  n <- max(vapply(v, length, integer(1)))
  lapply(v, rep_len, length.out = n)
}

#' The minor (rarer) allele base, given REF, ALT and the ALT frequency.
minor_allele_of <- function(ref, alt, af) {
  v  <- recycle_common(as.character(ref), as.character(alt), alt_is_minor(af))
  ifelse(is.na(v[[3]]), NA_character_, ifelse(v[[3]], v[[2]], v[[1]]))
}

#' The major (commoner) allele base.
major_allele_of <- function(ref, alt, af) {
  v  <- recycle_common(as.character(ref), as.character(alt), alt_is_minor(af))
  ifelse(is.na(v[[3]]), NA_character_, ifelse(v[[3]], v[[1]], v[[2]]))
}

#' Re-express an ALT-referenced slope per copy of the minor allele.
#'
#' Flips the sign wherever ALT is the major allele. NA where the reference
#' frequency is missing, because the minor allele is then unknown - reporting
#' the unflipped slope would silently mix conventions.
align_slope_to_minor <- function(slope, alt_minor) {
  v <- recycle_common(as.numeric(slope), alt_minor)
  ifelse(is.na(v[[2]]) | is.na(v[[1]]), NA_real_, v[[1]] * ifelse(v[[2]], 1, -1))
}

#' Minor-allele dosage from ALT dosage: `ploidy - alt_dosage` where ALT is the
#' major allele. Under trisomy 21 `ploidy` is 3, so the reflection is 3 - d.
minor_dosage <- function(alt_dosage, alt_minor, ploidy) {
  if (length(ploidy) != 1 || is.na(ploidy) || ploidy < 1)
    stop("ploidy must be a single number >= 1", call. = FALSE)
  bad <- !is.na(alt_dosage) & (alt_dosage < 0 | alt_dosage > ploidy)
  if (any(bad))
    stop("alt dosage outside [0, ", ploidy, "]: ", alt_dosage[which(bad)[1]], call. = FALSE)
  v <- recycle_common(as.numeric(alt_dosage), alt_minor)
  ifelse(is.na(v[[2]]) | is.na(v[[1]]), NA_real_,
         ifelse(v[[2]], v[[1]], ploidy - v[[1]]))
}

#' ALT allele frequency estimated from observed dosages at one variant.
#' Trisomic subjects carry three chr21 copies, so the denominator is `ploidy`.
alt_af_from_dosage <- function(alt_dosage, ploidy) {
  d <- alt_dosage[!is.na(alt_dosage)]
  if (!length(d)) return(NA_real_)
  mean(d) / ploidy
}

#' Minor-allele call for a set of variants, as a table.
#' @return data.table(maf, alt_is_minor, minor_allele, major_allele, maf_tie)
minor_call <- function(ref, alt, af) {
  data.table::data.table(
    maf          = maf_from_af(af),
    alt_is_minor = alt_is_minor(af),
    minor_allele = minor_allele_of(ref, alt, af),
    major_allele = major_allele_of(ref, alt, af),
    maf_tie      = af_is_tie(af))
}

#' Do two reference populations agree on which allele is minor?
#' NA when either frequency is missing, so "not checked" never reads as
#' "disagrees".
minor_allele_agrees <- function(af_a, af_b) {
  a <- alt_is_minor(af_a)
  b <- alt_is_minor(af_b)
  ifelse(is.na(a) | is.na(b), NA, a == b)
}
