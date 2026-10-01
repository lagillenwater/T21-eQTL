# 02_filter_genotypes.R
#
# Purpose: Build the genotype + cis-variant universe for the eQTL stage of
#          the pipeline.
#
#          Two changes from the previous version:
#          (a) variant source is GTEx v10 whole-blood ALLPAIRS (chr21
#              extract), not signif_pairs - so we get full cis-window
#              coverage per gene. We apply a nominal pval threshold to
#              keep the variant universe manageable.
#          (b) gene selection applies Hunter et al. (2023)'s own classification
#              rule BEFORE eQTL testing: a gene is pulled if it is eligible
#              (baseMean >= 30, Hunter et al.'s minimum coverage, not
#              repeat-flagged) AND
#              norm_padj < ALPHA_DE AND abs(norm_log2FC) >= DEVIATION_LFC_T2
#              (= log2(4/3); tier = 1 above DEVIATION_LFC = log2(1.5)), on
#              the ploidy-corrected scale. Genes whose
#              deviation is inside that threshold are not eQTL-tested - their
#              statistical significance is downstream of sample size, not
#              biological compensation, and asking the eQTL question for them
#              produces misleading detection calls.
#
#              The chr21-internal FDR-outlier test (OUTLIER_FDR, dev_z,
#              q_outlier) is retained as an ANNOTATION only. It does not
#              select genes and does not gate any lane. See the decision log
#              in docs/REPO_STATE.md and
#              docs/superpowers/plans/2026-08-31-tight-plan.md.
#
# Inputs:
#   - results/tables/deseq2_chr21_genes_both_analyses.csv
#   - data/processed/blacklisted_genes.csv
#   - data/GTEx_Analysis_v10_QTLs_GTEx_Analysis_v10_eQTL_all_associations_Whole_Blood.v10.allpairs.chr21.parquet
#   - data/processed/sample_metadata.csv
#   - data/chr21_ds_PASS.csv
#   - data/chr21_ctrl_PASS.csv
#
# Outputs:
#   - data/processed/eqtl_supported_genes.csv      (target gene roster)
#   - data/processed/eqtl_target_variants.csv      (cis variants per gene, with
#                                                   the minor-allele reference:
#                                                   GTEx af and gnomAD v4.1 AF)
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
source("scripts/lib/run.R"); run <- load_run(); th <- run$thresholds
source("scripts/lib/alleles.R")   # minor-allele referencing of the cis variants
source("scripts/lib/gnomad.R")    # gnomAD v4.1 ALT frequencies (cached tabix pulls)

cat(sprintf("=== T21-eQTL: Filter Genotypes for eQTL-Supported Genes [run: %s] ===\n\n", run$name))

# =============================================================================
# Constants - tunable filters
# =============================================================================

# Thresholds come from the run config (config/runs/<name>.R); the local names
# are kept so the rest of the script reads as before.
ALPHA_DE          <- th$alpha_de          # paper: padj < .01 after ploidy normalization
OUTLIER_FDR       <- th$outlier_fdr       # annotation only (dev_z, q_outlier)
DEVIATION_LFC     <- th$deviation_lfc     # tier 1: Hunter et al.'s FC >= 1.5 cut
DEVIATION_LFC_T2  <- th$deviation_lfc_t2  # tier 2: the cut target selection is made at
LOW_EXPR_BASEMEAN <- th$low_expr_basemean # Hunter et al.'s minimum coverage
GTEX_PVAL_KEEP    <- th$gtex_pval_keep    # nominal cis-eQTL pval cutoff in GTEx allpairs
source("scripts/lib/biotypes.R")   # TARGET_BIOTYPES: the chr21 target set (and
                                   # the chr21-internal null), the biotypes GTEx
                                   # whole blood tests

KNOWN_REPEAT_GENES <- c("RPS6KB1", "RPS27", "RPS27L", "RPS27P",
                        "IFNAR1", "IFNAR2", "TPTE", "BAGE", "DAB1")

# =============================================================================
# STEP 1: Identify target genes (DE + magnitude-filtered + paper filters)
# =============================================================================

cat("Step 1: Identifying target genes...\n")

chr21    <- fread(run$table("deseq2_chr21_genes_both_analyses.csv"))

blacklist_path <- "data/processed/blacklisted_genes.csv"
blacklist_genes <- if (file.exists(blacklist_path)) {
  fread(blacklist_path)$Gene_name
} else character(0)
high_repeat_genes <- unique(c(blacklist_genes, KNOWN_REPEAT_GENES))

n_chr21_before <- nrow(chr21)
chr21 <- chr21[Gene_type %in% TARGET_BIOTYPES]
cat(sprintf("  Restricted to %s: chr21 %d -> %d\n",
            paste(TARGET_BIOTYPES, collapse = " + "),
            n_chr21_before, nrow(chr21)))
print(chr21[, .N, by = Gene_type])

source("scripts/lib/chr21_threshold.R")

# Eligibility filters run BEFORE the null is estimated: log2FC variance scales
# with expression, so near-zero-count genes would otherwise set the scale. The
# cutoff is Hunter et al.'s absolute minimum, not a quantile of the target set,
# so adding lncRNA (mostly barely expressed in blood) cannot move it.
basemean_threshold <- LOW_EXPR_BASEMEAN
eligible <- chr21[baseMean >= basemean_threshold &
                    !(Gene_name %in% high_repeat_genes) &
                    !is.na(norm_log2FC)]
cat(sprintf("  Eligible for the null (expressed, non-repeat): %d of %d\n",
            nrow(eligible), nrow(chr21)))

null <- chr21_null(eligible$norm_log2FC)
cat(sprintf("  chr21 null: center %.4f  MAD %.4f  (n = %d)\n",
            null$center, null$scale, null$n))

eligible[, dev_z := robust_z(norm_log2FC, null)]
eligible[, q_outlier := outlier_fdr(dev_z)]
cat(sprintf("  Annotation only - FDR-outlier test at FDR < %.2f flags %d genes (effective k = %.2f)\n",
            OUTLIER_FDR, sum(eligible$q_outlier < OUTLIER_FDR, na.rm = TRUE),
            effective_k(eligible$dev_z, eligible$q_outlier, OUTLIER_FDR)))

# Composite rule: real deviation (padj) AND at least 1.5-fold (Hunter et al.'s
# effect-size cut on the ploidy-corrected scale).
# Two tiers, both eQTL-tested downstream. Tier 1 is Hunter's rule and is the
# primary result; tier 2 (>= 4/3-fold, < 1.5-fold) is the labelled secondary
# tier adopted 2026-09-01 so near-threshold genes are reported, not hidden.
target_genes <- eligible[!is.na(norm_padj) & norm_padj < ALPHA_DE &
                           abs(norm_log2FC) >= DEVIATION_LFC_T2]
if (nrow(target_genes) == 0) stop("no chr21 gene passes the deviation rule under run ", run$name, "; nothing to test", call. = FALSE)
target_genes[, tier := fifelse(abs(norm_log2FC) >= DEVIATION_LFC, 1L, 2L)]

target_genes[, gene_set := fifelse(norm_log2FC < 0,
                                   "DE_low_FC", "Sig_high_FC")]
target_genes[, ensembl_stable := sub("\\..*$", "", EnsemblID)]
target_genes[, observed_direction := sign(norm_log2FC)]
target_genes <- target_genes[, .(EnsemblID, Gene_name, raw_log2FC,
                                 norm_log2FC, norm_padj, gene_set,
                                 ensembl_stable, observed_direction, tier)]

n_low  <- sum(target_genes$gene_set == "DE_low_FC")
n_high <- sum(target_genes$gene_set == "Sig_high_FC")
cat(sprintf("  DE_low_FC genes passing magnitude filter:  %d\n", n_low))
cat(sprintf("  Sig_high_FC genes passing magnitude filter: %d\n", n_high))
cat(sprintf("  Total target genes for eQTL testing:        %d\n",
            nrow(target_genes)))

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

# ---- Positive-control genes (standalone control, tested in script 03) -------
# Controls v2. One positive control per tested deviating gene, drawn from the
# GTEx whole-blood eGenes (qval < positive_egene_qval) among the expressed,
# non-repeat chr21 genes that do NOT deviate, matched to that deviating gene
# on GTEx allelic fold change and expression level, one per locus (TSS at
# least positive_min_separation apart), and never a gene in a GRCh38 false
# duplication (data/grch38_false_duplication_genes_chr21.csv), where read
# mapping and genotyping are unreliable. Same variant pull here and the same
# within-T21 permutation test in script 03 as the deviating genes, under
# gene_set == "positive_control"; scripts 03 and 04 keep them out of every
# main-result table. They show what the test does on real eQTLs of the size
# it is asked to find, not on the strongest eQTLs on the chromosome.
source("scripts/lib/eqtl_controls.R")
EGENES_PATH  <- "data/Whole_Blood.v10.eGenes.txt.gz"
FALSE_DUP    <- fread("data/grch38_false_duplication_genes_chr21.csv")
egenes <- fread(EGENES_PATH)[gene_chr == "chr21"]
egenes[, `:=`(ensembl_stable = sub("\\..*$", "", gene_id),
              tss_egenes = fifelse(strand == "+", as.numeric(gene_start), as.numeric(gene_end)))]
gtex_min_p <- gtex[startsWith(variant_id, "chr21_") & !is.na(pval_nominal),
                   .(gtex_min_p = min(pval_nominal)), by = ensembl_stable]
eligible[, ensembl_stable := sub("\\..*$", "", EnsemblID)]
# Targets: the deviating genes that will be tested (at least one GTEx cis
# variant at the pval cut). Candidates also keep positive_dev_separation
# (100 kb) from every deviating gene's TSS, so a control is never an
# overlapping or antisense partner of a deviating gene.
testable <- gtex_min_p[gtex_min_p <= GTEX_PVAL_KEEP, ensembl_stable]
# eGene strength: the deviating genes are all strong GTEx eGenes, so a
# candidate must be one too (positive_egene_qval) and carry at least
# positive_min_variants variants at the pval cut, or its "effect" is one
# GTEx barely established and a one-variant set in T21.
n_at_cut <- gtex[startsWith(variant_id, "chr21_") & !is.na(pval_nominal) & pval_nominal <= GTEX_PVAL_KEEP,
                 .(n_at_cut = .N), by = ensembl_stable]
pc_targets <- merge(eligible[ensembl_stable %in% intersect(target_genes$ensembl_stable, testable),
                             .(Gene_name, ensembl_stable, baseMean)],
                    egenes[, .(ensembl_stable, abs_afc = abs(afc))], by = "ensembl_stable")
dev_tss <- egenes[ensembl_stable %in% target_genes$ensembl_stable, tss_egenes]
pc_candidates <- merge(eligible[!ensembl_stable %in% target_genes$ensembl_stable &
                                  !Gene_name %in% FALSE_DUP$gene,
                                .(ensembl_stable, EnsemblID, Gene_name, baseMean, raw_log2FC, norm_log2FC, norm_padj)],
                       egenes[qval < th$positive_egene_qval,
                              .(ensembl_stable, abs_afc = abs(afc), tss = tss_egenes, gtex_qval = qval)],
                       by = "ensembl_stable")
pc_candidates <- merge(pc_candidates, n_at_cut, by = "ensembl_stable")[n_at_cut >= th$positive_min_variants]
near_dev <- apply(abs(outer(pc_candidates$tss, dev_tss, "-")) < th$positive_dev_separation, 1, any)
pc_candidates <- pc_candidates[!near_dev]
cat(sprintf("  Positive-control pool: %d eGenes (q < %.2f) among %d eligible non-deviating genes; %d false-duplication genes excluded by name\n",
            nrow(pc_candidates), th$positive_egene_qval, sum(!eligible$ensembl_stable %in% target_genes$ensembl_stable),
            sum(eligible$Gene_name %in% FALSE_DUP$gene)))
matching <- match_positive_controls(pc_targets, pc_candidates, th$positive_min_separation)
fwrite(matching, run$table("positive_control_matching.csv"))
positive_controls <- merge(pc_candidates[ensembl_stable %in% matching$ensembl_stable,
                                         .(ensembl_stable, EnsemblID, Gene_name, raw_log2FC, norm_log2FC, norm_padj)],
                           gtex_min_p, by = "ensembl_stable", all.x = TRUE)
positive_controls <- merge(positive_controls, matching[, .(ensembl_stable, matched_to = target_gene)], by = "ensembl_stable")
positive_controls[, `:=`(gene_set = "positive_control",
                         observed_direction = sign(norm_log2FC),
                         tier = NA_integer_)]
cat(sprintf("  Positive-control genes (matched GTEx eGenes, one per locus, non-deviating): %d for %d targets\n",
            nrow(positive_controls), nrow(pc_targets)))
print(matching[, .(target_gene, target_abs_afc = round(target_abs_afc, 2), target_baseMean = round(target_baseMean),
                   control = Gene_name, abs_afc = round(abs_afc, 2), baseMean = round(baseMean),
                   tss_mb = round(tss / 1e6, 2), match_distance = round(match_distance, 2))])
target_genes <- rbind(target_genes, positive_controls, fill = TRUE)

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

# TSS per gene (GTEx: tss_distance = POS - TSS), carried on the roster so
# script 03 can pair each deviating gene with a distant decoy variant set and
# script 11 can place it on the chr21 map. Taken from every allpairs row of
# the gene, not only those passing the p cut, so a gene with GTEx rows but no
# variant at the cut still has a position.
tss_rows <- gtex[ensembl_stable %in% target_genes$ensembl_stable &
                   startsWith(variant_id, "chr21_")]
tss_rows[, POS := as.integer(tstrsplit(variant_id, "_", fixed = TRUE)[[2]])]
target_genes <- merge(target_genes,
                      gene_tss(tss_rows[, .(ensembl_stable, POS, tss_distance)]),
                      by = "ensembl_stable", all.x = TRUE)
fwrite(target_genes, run$processed("eqtl_supported_genes.csv"))

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

# ---- Minor-allele reference for every cis variant --------------------------
# Which of REF/ALT is the rarer allele in the typical population. Genotype
# stays ALT-coded downstream; these columns are what lets scripts 03, 04 and
# 11 report a direction per MINOR allele, the allele eQTL effects are normally
# reported against (scripts/lib/alleles.R). GTEx's own `af` is primary - it is
# the population the cis-eQTL was called in - and gnomAD v4.1 genomes is the
# independent check on the same call.

target_variants[, c("gtex_maf", "alt_is_minor", "minor_allele",
                    "major_allele", "maf_tie") :=
                  minor_call(REF, ALT, af)]
stopifnot(!anyNA(target_variants$alt_is_minor))   # every allpairs row carries af

gnomad <- ensure_gnomad_af(target_variants$POS)
gnomad_status <- attr(gnomad, "gnomad_status")
if (!identical(gnomad_status, "ok"))
  warning("gnomAD frequencies incomplete (", gnomad_status,
          "); the gnomad_* columns are NA where they could not be fetched. ",
          "Re-run scripts/fetch_gnomad_af.R once the network is available.",
          call. = FALSE)

if (nrow(gnomad)) {
  setnames(gnomad, c("AF", "AF_nfe", "filter"),
           c("gnomad_af", "gnomad_af_nfe", "gnomad_filter"))
  target_variants <- merge(target_variants, gnomad, by = c("POS", "REF", "ALT"), all.x = TRUE)
} else {
  target_variants[, c("gnomad_af", "gnomad_af_nfe", "gnomad_filter") :=
                    .(NA_real_, NA_real_, NA_character_)]
}
target_variants[, `:=`(
  gnomad_maf            = maf_from_af(gnomad_af),
  gnomad_alt_is_minor   = alt_is_minor(gnomad_af),
  minor_concordant      = minor_allele_agrees(af, gnomad_af),
  minor_concordant_nfe  = minor_allele_agrees(af, gnomad_af_nfe))]

n_var <- nrow(target_variants)
cat(sprintf("  ALT is the minor allele in GTEx whole blood: %d / %d (%.1f%%)\n",
            sum(target_variants$alt_is_minor), n_var,
            100 * mean(target_variants$alt_is_minor)))
cat(sprintf("  gnomAD v4.1 frequency found:                %d / %d (%.1f%%)\n",
            sum(!is.na(target_variants$gnomad_af)), n_var,
            100 * mean(!is.na(target_variants$gnomad_af))))
cat(sprintf("  GTEx and gnomAD agree on the minor allele:  %d / %d matched (%.1f%%)\n",
            sum(target_variants$minor_concordant, na.rm = TRUE),
            sum(!is.na(target_variants$minor_concordant)),
            100 * mean(target_variants$minor_concordant, na.rm = TRUE)))

fwrite(target_variants, run$processed("eqtl_target_variants.csv"))

target_pos <- sort(unique(target_variants$POS))
cat(sprintf("  Unique positions to scan in genotype CSVs: %d\n",
            length(target_pos)))

# =============================================================================
# STEP 3: DESeq2 sample roster (subject base ID -> karyotype)
# =============================================================================

cat("\nStep 3: Loading DESeq2 sample roster...\n")

meta <- run$cohort()
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

geno_out <- run$processed("genotypes_filtered.csv")
pos_sig  <- run$processed("genotype_target_positions.txt")
# The cache key covers the target positions AND the cohort subjects, because
# the stream also column-filters to the run cohort.
cache_key <- c(as.character(target_pos), "--subjects--", sort(unique(meta$subject_id)))
if (file.exists(geno_out) && file.exists(pos_sig) &&
    identical(readLines(pos_sig), cache_key)) {
  cat("  Genotype extract for these target positions already exists; skipping the stream.\n")
  out <- fread(geno_out)
} else {

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
fwrite(out, geno_out)
writeLines(cache_key, pos_sig)
cat(sprintf("\n  Wrote %s (%d rows)\n", geno_out, nrow(out)))
}

# =============================================================================
# STEP 6: Verification summary
# =============================================================================

cat("\n=== Verification ===\n")
n_low_fc  <- uniqueN(
  target_genes$ensembl_stable[target_genes$gene_set == "DE_low_FC"])
n_high_fc <- uniqueN(
  target_genes$ensembl_stable[target_genes$gene_set == "Sig_high_FC"])
cat(sprintf("Target genes: %d  (DE_low_FC: %d, Sig_high_FC: %d, positive_control: %d)\n",
            uniqueN(target_genes$ensembl_stable), n_low_fc, n_high_fc,
            sum(target_genes$gene_set == "positive_control")))
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

writeLines(capture.output(sessionInfo()), run$processed("genotype_filter_session_info.txt"))

stopifnot(
  file.exists(run$processed("genotype_filter_session_info.txt")),
  file.exists(run$processed("eqtl_supported_genes.csv")),
  file.exists(run$processed("eqtl_target_variants.csv")),
  file.exists(geno_out)
)

cat("\n=== Filter complete ===\n")

# =============================================================================
# CHANGELOG
# =============================================================================
# 2026-08-31  REPLACED the cohort-SD magnitude filter
#             (abs(norm_log2FC) >= MAGNITUDE_THRESHOLD * sd(non-chr21 log2FC),
#             threshold 1.0) with an FDR-controlled robust outlier test against
#             a chr21-internal median/MAD null (scripts/lib/chr21_threshold.R,
#             OUTLIER_FDR = 0.10).
#             Reason: ploidy normalization does not act on diploid genes
#             (mean |raw - norm| 0.0048 off chr21 vs 0.583 on it), so their
#             spread measured a different quantity; and a 1-SD cut selects the
#             top ~third of any distribution (18.3% of non-chr21 genes cleared
#             it themselves, vs 20.6% of chr21 - binomial p = 0.26). The null is
#             now estimated AFTER the expression and repeat filters, because
#             log2FC variance scales with counts.
#             Spec: docs/METHODS_SPEC_threshold_and_eqtl_controls.md
#
# 2026-08-31  REPLACED the FDR-outlier test in the target-gene rule with
#             Hunter et al.'s own classification: deviating = padj < 0.01 AND
#             abs(norm_log2FC) >= log2(1.5) (DEVIATION_LFC). dev_z / q_outlier
#             are retained as annotation columns on the eligible table but no
#             longer drive target-gene selection.
#             Spec: docs/superpowers/plans/2026-08-31-tight-plan.md (Task A)
#
# 2026-08-31  REWROTE the header and the OUTLIER_FDR comment, which still
#             described the FDR-outlier test as the gene-selection rule after
#             the rule itself had been replaced by Hunter's padj + 1.5-fold
#             cut. No behaviour change: the code already selected on
#             ALPHA_DE + DEVIATION_LFC. The FDR-outlier print is now labelled
#             annotation-only so the two cannot be confused again.
# 2026-09-01  ADDED tier 2 (DEVIATION_LFC_T2 = log2(4/3)): the target set now
#             carries both tiers with a `tier` column; tier 1 (log2(1.5),
#             Hunter) remains the primary result. Adopted after review of the
#             volcano figure showed genes with crushing padj within 0.02-0.09
#             of the tier-1 line; reported as a labelled secondary tier rather
#             than moving the pre-registered primary threshold.
# 2026-09-04  WIDENED the chr21 target set from protein-coding only to
#             protein-coding, lncRNA and pseudogene biotypes (TARGET_BIOTYPES in
#             scripts/lib/biotypes.R), the biotypes GTEx whole blood tests. REPLACED the q20 baseMean
#             low-expression cutoff with Hunter et al.'s absolute minimum
#             (LOW_EXPR_BASEMEAN = 30): a quantile of the target set would have
#             dropped from 25.1 to 4.1 once lncRNA joined it, admitting genes
#             with a handful of counts. On the protein-coding set the absolute
#             cutoff flags 34 genes against the quantile's 32 (a superset).
# 2026-09-10  ADDED positive-control genes to the roster (gene_set ==
#             "positive_control"): the n_positive_controls strongest GTEx
#             whole-blood eGenes among expressed, non-repeat, non-deviating
#             chr21 genes (scripts/lib/eqtl_controls.R), pulled through the
#             same variant and genotype extraction. ADDED gtex_min_p and tss
#             columns to eqtl_supported_genes.csv; the roster is now written
#             after the GTEx pull. Deviating-gene selection is unchanged.
