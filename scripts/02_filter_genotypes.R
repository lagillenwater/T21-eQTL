# 09_filter_genotypes.R
#
# Purpose: Build the genotype + cis-variant universe for the eQTL stage of
#          the pipeline.
#
#          Two changes from the previous version:
#          (a) variant source is GTEx v10 whole-blood ALLPAIRS (chr21
#              extract), not signif_pairs - so we get full cis-window
#              coverage per gene. We apply a nominal pval threshold to
#              keep the variant universe manageable.
#          (b) gene selection now applies a within-cohort SD magnitude
#              filter BEFORE eQTL testing. Only chr21 DE genes whose
#              ploidy-corrected deviation exceeds MAGNITUDE_THRESHOLD
#              cohort SDs are pulled. Genes with deviations within
#              typical cohort variation are not eQTL-tested - their
#              statistical significance is downstream of sample size,
#              not biological compensation, and asking the eQTL question
#              for them produces misleading "explained" calls.
#
# Inputs:
#   - results/tables/deseq2_chr21_combined.csv
#   - results/tables/deseq2_all_genes_ploidy_normalized.csv
#   - data/processed/blacklisted_genes.csv
#   - data/GTEx_Analysis_v10_QTLs_GTEx_Analysis_v10_eQTL_all_associations_Whole_Blood.v10.allpairs.chr21.parquet
#   - data/processed/sample_metadata.csv
#   - data/chr21_ds_PASS.csv
#   - data/chr21_ctrl_PASS.csv
#
# Outputs:
#   - data/processed/eqtl_supported_genes.csv      (target gene roster)
#   - data/processed/eqtl_target_variants.csv      (cis variants per gene)
#   - data/processed/genotypes_filtered.csv        (HTP genotypes at those
#                                                   variants, T21+Control)
#   - data/processed/genotype_filter_session_info.txt
#
# Date: 2026-05-04

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(arrow)
})

set.seed(42)

cat("=== T21-eQTL: Filter Genotypes for eQTL-Supported Genes ===\n\n")

# =============================================================================
# Constants - tunable filters
# =============================================================================

ALPHA_DE          <- 0.01    # paper: padj < .01 after ploidy normalization
MAGNITUDE_THRESHOLD <- 1.0   # |norm_log2FC| / cohort_sd >= this to be
                             # eQTL-tested. 1.0 = above typical cohort
                             # variation. Raise to 1.5 or 2.0 for stricter.
LOW_EXPR_QUANT    <- 0.20    # paper's 2nd-quintile baseMean filter
GTEX_PVAL_KEEP    <- 1e-4    # nominal cis-eQTL pval cutoff in GTEx allpairs
                             # (~ matches the effective signif_pairs cutoff)
RESTRICT_TO_PROTEIN_CODING <- TRUE   # restrict both target chr21 set AND
                                     # the non-chr21 cohort-noise reference
                                     # to protein-coding genes

KNOWN_REPEAT_GENES <- c("RPS6KB1", "RPS27", "RPS27L", "RPS27P",
                        "IFNAR1", "IFNAR2", "TPTE", "BAGE", "DAB1")

# =============================================================================
# STEP 1: Identify target genes (DE + magnitude-filtered + paper filters)
# =============================================================================

cat("Step 1: Identifying target genes...\n")

chr21    <- fread("results/tables/deseq2_chr21_combined.csv")
all_lfc  <- fread(
  "results/tables/deseq2_all_genes_ploidy_normalized.csv")

blacklist_path <- "data/processed/blacklisted_genes.csv"
blacklist_genes <- if (file.exists(blacklist_path)) {
  fread(blacklist_path)$Gene_name
} else character(0)
high_repeat_genes <- unique(c(blacklist_genes, KNOWN_REPEAT_GENES))

if (RESTRICT_TO_PROTEIN_CODING) {
  n_chr21_before <- nrow(chr21)
  chr21   <- chr21[Gene_type == "protein_coding"]
  all_lfc <- all_lfc[Gene_type == "protein_coding"]
  cat(sprintf("  Restricted to protein-coding: chr21 %d -> %d\n",
              n_chr21_before, nrow(chr21)))
}

# Cohort-noise SD from non-chr21 ploidy-normalized log2FC (the empirical
# null distribution for "deviation from expected ploidy"). Restricted to
# protein-coding too if RESTRICT_TO_PROTEIN_CODING is TRUE.
cohort_sd <- sd(
  all_lfc[Chr != "chr21" & !is.na(log2FoldChange), log2FoldChange])
cat(sprintf("  Cohort-noise SD (non-chr21%s): %.3f\n",
            if (RESTRICT_TO_PROTEIN_CODING) ", protein-coding" else "",
            cohort_sd))
cat(sprintf("  Magnitude threshold: %.1f cohort SDs (= %.3f log2FC)\n",
            MAGNITUDE_THRESHOLD, MAGNITUDE_THRESHOLD * cohort_sd))

basemean_threshold <- quantile(chr21$baseMean, LOW_EXPR_QUANT, na.rm = TRUE)

target_genes <- chr21[
  !is.na(norm_padj) & norm_padj < ALPHA_DE &
    !is.na(norm_log2FC) &
    abs(norm_log2FC) >= MAGNITUDE_THRESHOLD * cohort_sd &
    baseMean >= basemean_threshold &
    !(Gene_name %in% high_repeat_genes)]

target_genes[, gene_set := fifelse(norm_log2FC < 0,
                                   "DE_low_FC", "Sig_high_FC")]
target_genes[, ensembl_stable := sub("\\..*$", "", EnsemblID)]
target_genes[, observed_direction := sign(norm_log2FC)]
target_genes <- target_genes[, .(EnsemblID, Gene_name, raw_log2FC,
                                 norm_log2FC, norm_padj, gene_set,
                                 ensembl_stable, observed_direction)]

n_low  <- sum(target_genes$gene_set == "DE_low_FC")
n_high <- sum(target_genes$gene_set == "Sig_high_FC")
cat(sprintf("  DE_low_FC genes passing magnitude filter:  %d\n", n_low))
cat(sprintf("  Sig_high_FC genes passing magnitude filter: %d\n", n_high))
cat(sprintf("  Total target genes for eQTL testing:        %d\n",
            nrow(target_genes)))

fwrite(target_genes, "data/processed/eqtl_supported_genes.csv")

# =============================================================================
# STEP 2: Pull GTEx whole-blood ALLPAIRS cis variants for the target genes
# =============================================================================

cat("\nStep 2: Loading GTEx whole-blood allpairs (chr21)...\n")

allpairs_path <- paste0("data/GTEx_Analysis_v10_QTLs_GTEx_Analysis_v10_eQTL",
                        "_all_associations_Whole_Blood.v10.allpairs.chr21",
                        ".parquet")
gtex <- as.data.table(read_parquet(allpairs_path))
cat(sprintf("  allpairs rows loaded: %d (genes: %d, variants: %d)\n",
            nrow(gtex), uniqueN(gtex$gene_id), uniqueN(gtex$variant_id)))

gtex[, ensembl_stable := sub("\\..*$", "", gene_id)]

target_variants <- gtex[ensembl_stable %in% target_genes$ensembl_stable &
                        startsWith(variant_id, "chr21_") &
                        !is.na(pval_nominal) &
                        pval_nominal <= GTEX_PVAL_KEEP]
cat(sprintf("  After pval_nominal <= %.0e filter: %d rows\n",
            GTEX_PVAL_KEEP, nrow(target_variants)))

# Parse variant_id: chr21_POS_REF_ALT_b38
parsed <- tstrsplit(target_variants$variant_id, "_", fixed = TRUE)
target_variants[, `:=`(
  CHROM = parsed[[1]],
  POS   = as.integer(parsed[[2]]),
  REF   = parsed[[3]],
  ALT   = parsed[[4]]
)]

# Attach gene metadata (some variants may map to multiple genes - keep all)
target_variants <- merge(
  target_variants,
  target_genes[, .(ensembl_stable, Gene_name, gene_set,
                   raw_log2FC, norm_log2FC, observed_direction)],
  by = "ensembl_stable", allow.cartesian = TRUE
)

cat(sprintf("  cis variants for target genes: %d\n", nrow(target_variants)))
cat(sprintf("  Unique target variants (POS,REF,ALT): %d\n",
            uniqueN(target_variants[, .(POS, REF, ALT)])))
cat(sprintf("  Target genes covered by GTEx allpairs: %d / %d\n",
            uniqueN(target_variants$ensembl_stable),
            uniqueN(target_genes$ensembl_stable)))

fwrite(target_variants, "data/processed/eqtl_target_variants.csv")

target_pos <- sort(unique(target_variants$POS))
cat(sprintf("  Unique positions to scan in genotype CSVs: %d\n",
            length(target_pos)))

# =============================================================================
# STEP 3: DESeq2 sample roster (subject base ID -> karyotype)
# =============================================================================

cat("\nStep 3: Loading DESeq2 sample roster...\n")

meta <- fread("data/processed/sample_metadata.csv")
meta[, subject_id := sub("[A-Z][0-9]*$", "", LabID)]   # strip A/A2/B/B2/...

stopifnot(all(meta$Karyotype %in% c("T21", "Control")))
cat(sprintf("  Samples in DESeq2 cohort: %d  (T21: %d, Control: %d)\n",
            nrow(meta), sum(meta$Karyotype == "T21"),
            sum(meta$Karyotype == "Control")))

deseq_subjects <- unique(meta$subject_id)
subject_karyo <- setNames(meta$Karyotype[!duplicated(meta$subject_id)],
                          meta$subject_id[!duplicated(meta$subject_id)])

# =============================================================================
# STEP 4: Stream-filter genotype CSVs to target positions and DESeq2 samples
# =============================================================================

cat("\nStep 4: Filtering genotype CSVs (stream via awk)...\n")

# Write target positions to a temp file for awk to load into a hash
pos_tmp <- tempfile(fileext = ".txt")
writeLines(as.character(target_pos), pos_tmp)

read_filtered_genotypes <- function(geno_path, pos_file, deseq_subjects,
                                    subject_karyo, expected_karyo) {
  cat(sprintf("  -> %s\n", geno_path))

  awk_cmd <- sprintf(
    paste0("awk -F',' 'BEGIN{while((getline l < \"%s\")>0) p[l]=1} ",
           "NR==1 || ($2 in p)' %s"),
    pos_file, geno_path
  )

  dt <- fread(cmd = awk_cmd, sep = ",", header = TRUE,
              showProgress = FALSE)

  cat(sprintf("     rows after POS filter: %d\n", nrow(dt)))

  meta_cols <- c("CHROM", "POS", "ID", "REF", "ALT",
                 "QUAL", "FILTER", "INFO", "FORMAT")
  geno_cols <- setdiff(colnames(dt), meta_cols)

  sample_subjects <- sub("[A-Z][0-9]*$", "", geno_cols)
  keep_mask <- sample_subjects %in% deseq_subjects &
    subject_karyo[sample_subjects] == expected_karyo
  keep_cols <- geno_cols[keep_mask]

  cat(sprintf("     samples kept (matched %s in DESeq2): %d\n",
              expected_karyo, length(keep_cols)))

  if (length(keep_cols) == 0) return(NULL)

  dt <- dt[, c(meta_cols, keep_cols), with = FALSE]
  long <- melt(dt,
               id.vars = meta_cols,
               measure.vars = keep_cols,
               variable.name = "lab_id",
               value.name = "geno_field",
               variable.factor = FALSE)

  long[, `:=`(
    subject_id = sub("[A-Z][0-9]*$", "", lab_id),
    karyotype  = expected_karyo
  )]

  long
}

ds_long   <- read_filtered_genotypes("data/chr21_ds_PASS.csv", pos_tmp,
                                     deseq_subjects, subject_karyo, "T21")
ctrl_long <- read_filtered_genotypes("data/chr21_ctrl_PASS.csv", pos_tmp,
                                     deseq_subjects, subject_karyo, "Control")

geno_long <- rbindlist(list(ds_long, ctrl_long), use.names = TRUE)

cat(sprintf("\n  Combined genotype rows (variant x sample): %d\n",
            nrow(geno_long)))

# =============================================================================
# STEP 5: Restrict to (POS, REF, ALT) matches and parse GT to alt dosage
# =============================================================================

cat("\nStep 5: Matching REF/ALT to GTEx variants and parsing dosage...\n")

# Drop multi-allelic ALT rows (we keep biallelic exact matches only for the
# directional test - GTEx variant_id encodes a single ALT)
geno_long <- geno_long[!grepl(",", ALT)]

# Inner-join to target variants on (POS, REF, ALT) - this also drops positions
# where GTEx variant_id ALT does not match the HTP call
geno_long <- merge(
  geno_long,
  unique(target_variants[, .(variant_id, POS, REF, ALT)]),
  by = c("POS", "REF", "ALT"), all = FALSE
)

cat(sprintf("  Rows after REF/ALT match: %d\n", nrow(geno_long)))

# Parse GT (first colon-subfield); count number of alt alleles ('1') in
# phased ('|') or unphased ('/') genotype strings. Trisomic chr21 in T21
# samples can yield 3 alleles (e.g., 0|0|1 -> 1, 1|1|1 -> 3).
parse_alt_dosage <- function(geno_field) {
  gt <- sub(":.*$", "", geno_field)
  alleles <- strsplit(gt, "[|/]", fixed = FALSE)
  vapply(alleles, function(a) {
    a <- a[a != "."]
    if (length(a) == 0L) NA_integer_ else sum(a == "1")
  }, integer(1))
}

geno_long[, alt_dosage := parse_alt_dosage(geno_field)]

n_missing <- sum(is.na(geno_long$alt_dosage))
cat(sprintf("  Missing/unparseable genotypes: %d (%.2f%%)\n",
            n_missing, 100 * n_missing / nrow(geno_long)))

# Final tidy table
out <- geno_long[, .(
  variant_id, CHROM, POS, REF, ALT,
  lab_id, subject_id, karyotype, alt_dosage
)]

setorder(out, POS, karyotype, subject_id)
fwrite(out, "data/processed/genotypes_filtered.csv")

cat(sprintf("\n  Wrote data/processed/genotypes_filtered.csv (%d rows)\n",
            nrow(out)))

# =============================================================================
# STEP 6: Verification summary
# =============================================================================

cat("\n=== Verification ===\n")
n_low_fc  <- uniqueN(
  target_genes$ensembl_stable[target_genes$gene_set == "DE_low_FC"])
n_high_fc <- uniqueN(
  target_genes$ensembl_stable[target_genes$gene_set == "Sig_high_FC"])
cat(sprintf("Target genes: %d  (DE_low_FC: %d, Sig_high_FC: %d)\n",
            uniqueN(target_genes$ensembl_stable), n_low_fc, n_high_fc))
cat(sprintf("Variants tested in HTP: %d / %d GTEx targets\n",
            uniqueN(out$variant_id), uniqueN(target_variants$variant_id)))
cat(sprintf("Subjects with genotypes: %d  (T21: %d, Control: %d)\n",
            uniqueN(out$subject_id),
            uniqueN(out$subject_id[out$karyotype == "T21"]),
            uniqueN(out$subject_id[out$karyotype == "Control"])))

dosage_range <- out[, .(min = min(alt_dosage, na.rm = TRUE),
                        max = max(alt_dosage, na.rm = TRUE)),
                    by = karyotype]
print(dosage_range)
cat("(Expect Control max <= 2; T21 max <= 3 on chr21)\n")

writeLines(capture.output(sessionInfo()),
           "data/processed/genotype_filter_session_info.txt")

cat("\n=== Filter complete ===\n")
