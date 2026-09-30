# verify_hunter_overlap.R
#
# Purpose: Check this analysis's deviating chr21 genes against Hunter et al.
#          (2023)'s own deposited DESeq2 results, in both directions, so the
#          manuscript's consistency claims can be checked rather than asserted.
#
# What it reproduces:
#   1. Hunter's lower-than-expected chr21 set, by re-applying their criterion
#      (chr21, padj < 0.01, log2FC < 0 on the ploidy-modified table) to their
#      deposited results. The paper names only the three that eQTLs explained
#      (CLIC6, ITSN1, C2CD2); the rest of the set is recovered here, so treat
#      the unnamed genes as reconstructed rather than cited.
#   2. Their higher-than-expected chr21 set, which the paper does not discuss.
#   3. Where each of their genes lands in this run's lane table, and whether
#      each of this run's deviating genes was flagged by them - including
#      whether it was expressed in their LCLs at all.
#
# Multi-mapping: three p-arm genes share one identical (baseMean, log2FC,
# padj) triple across several loci, the signature of one multi-mapped block
# counted repeatedly. Any statistics triple shared by more than one gene
# symbol is dropped here, which is what leaves Hunter's five.
#
# Cohorts differ (one LCL family trio vs 302 T21 whole-blood subjects), so
# neither overlap nor non-overlap is a reproduction test. This script reports
# what the two gene sets do, nothing more.
#
# Inputs:  data/hunter/rna_modified_deseq2.results.csv    (downloaded if absent)
#          data/hunter/refseq_to_common_id.txt            (downloaded if absent)
#          results/runs/<run>/tables/chr21_lane_assignments.csv
# Outputs: results/runs/<run>/tables/hunter_comparison.csv
#          results/runs/<run>/tables/hunter_comparison_reverse.csv
#
# Usage:   Rscript scripts/verify_hunter_overlap.R --run adjusted

suppressPackageStartupMessages(library(data.table))
source("scripts/lib/run.R"); run <- load_run()

ALPHA_HUNTER <- 0.01   # Hunter et al.'s padj threshold, the same cut as ours

# Deposited with the paper: Dowell-Lab/DS_Normalization, tag BMCbio2023.
HUNTER_BASE  <- "https://raw.githubusercontent.com/Dowell-Lab/DS_Normalization/BMCbio2023"
HUNTER_DIR   <- file.path("data", "hunter")
HUNTER_FILES <- c(rna_modified_deseq2.results.csv = "results_table/rna_modified_deseq2.results.csv",
                  refseq_to_common_id.txt         = "annotation/refseq_to_common_id.txt")

cat(sprintf("=== Hunter et al. (2023) comparison [run: %s] ===\n\n", run$name))

# =============================================================================
# STEP 1: Hunter's deposited results
# =============================================================================

dir.create(HUNTER_DIR, recursive = TRUE, showWarnings = FALSE)
for (nm in names(HUNTER_FILES)) {
  dest <- file.path(HUNTER_DIR, nm)
  if (!file.exists(dest)) {
    cat("  downloading", nm, "...\n")
    ok <- try(download.file(file.path(HUNTER_BASE, HUNTER_FILES[[nm]]), dest,
                            mode = "wb", quiet = TRUE), silent = TRUE)
    if (inherits(ok, "try-error")) {
      unlink(dest)
      stop("could not download ", HUNTER_FILES[[nm]], " - fetch it by hand into ",
           HUNTER_DIR, call. = FALSE)
    }
  }
}

hunter <- fread(file.path(HUNTER_DIR, "rna_modified_deseq2.results.csv"))
refseq <- fread(file.path(HUNTER_DIR, "refseq_to_common_id.txt"), header = FALSE,
                col.names = c("refseq", "Gene_name"))
hunter[, refseq := Row.names]
hunter <- merge(hunter, refseq, by = "refseq", all.x = TRUE)
cat(sprintf("Step 1: Hunter rows: %d (chr21: %d)\n", nrow(hunter), sum(hunter$chr == "chr21")))

sig <- unique(hunter[chr == "chr21" & !is.na(padj) & padj < ALPHA_HUNTER,
                     .(refseq, Gene_name, baseMean, log2FoldChange, padj, start, end)])
cat(sprintf("  chr21, padj < %.2f, distinct loci: %d\n", ALPHA_HUNTER, nrow(sig)))

# One multi-mapped block counted at several loci shows up as one statistics
# triple carrying several gene symbols.
sig[, n_symbols := uniqueN(Gene_name), by = .(baseMean, log2FoldChange, padj)]
multi <- sig[n_symbols > 1]
if (nrow(multi))
  cat(sprintf("  dropped as multi-mapped (identical statistics at %d loci): %s\n",
              nrow(multi), paste(sort(unique(multi$Gene_name)), collapse = ", ")))
sig <- sig[n_symbols == 1]

sig[, hunter_direction := fifelse(log2FoldChange < 0, "lower", "higher")]
setorder(sig, hunter_direction, log2FoldChange)
cat(sprintf("  Hunter deviating chr21 genes: %d lower, %d higher than expected\n",
            sum(sig$hunter_direction == "lower"), sum(sig$hunter_direction == "higher")))
print(sig[, .(Gene_name, hunter_direction, log2FC = round(log2FoldChange, 2),
              padj = signif(padj, 2))])

# The paper names these three of the lower set; recovering them is the check
# that the reconstruction of the rest is right.
NAMED_IN_PAPER <- c("CLIC6", "ITSN1", "C2CD2")
recovered <- NAMED_IN_PAPER %in% sig[hunter_direction == "lower", Gene_name]
cat(sprintf("\n  Genes the paper names as lower-than-expected and eQTL-explained: %s\n",
            paste(sprintf("%s %s", NAMED_IN_PAPER, ifelse(recovered, "OK", "MISSING")),
                  collapse = ", ")))
if (!all(recovered))
  warning("the filter did not recover every gene the paper names; the ",
          "reconstruction of the unnamed genes cannot be trusted", call. = FALSE)

# =============================================================================
# STEP 2: This run's chr21 lane table
# =============================================================================

lanes <- fread(run$table("chr21_lane_assignments.csv"))
ours  <- lanes[sig_lane %in% c("DE_high", "DE_low")]
ours[, our_direction := fifelse(sig_lane == "DE_high", "higher", "lower")]
cat(sprintf("\nStep 2: this run's deviating chr21 genes: %d (%d higher, %d lower)\n",
            nrow(ours), sum(ours$our_direction == "higher"),
            sum(ours$our_direction == "lower")))

# =============================================================================
# STEP 3: Hunter's genes -> our classification
# =============================================================================

cmp <- merge(sig[, .(Gene_name, hunter_direction, hunter_log2FC = log2FoldChange,
                     hunter_padj = padj)],
             lanes[, .(Gene_name, our_baseMean = baseMean, our_norm_log2FC = norm_log2FC,
                       our_norm_padj = norm_padj, our_sig_lane = sig_lane,
                       our_tier = tier, our_eqtl_lane = eqtl_lane)],
             by = "Gene_name", all.x = TRUE)
cmp[, our_direction := fcase(is.na(our_sig_lane), NA_character_,
                             our_sig_lane == "DE_high", "higher",
                             our_sig_lane == "DE_low",  "lower",
                             default = NA_character_)]
cmp[, verdict := fcase(
  is.na(our_sig_lane),                       "not in our chr21 target set",
  is.na(our_direction),                      paste0("not deviating here (", our_sig_lane, ")"),
  our_direction == hunter_direction,         "deviating in both, same direction",
  default =                                  "deviating in both, OPPOSITE direction")]
setorder(cmp, hunter_direction, hunter_log2FC)
fwrite(cmp, run$table("hunter_comparison.csv"))

cat("\nStep 3: where Hunter's deviating genes land here\n")
print(cmp[, .(Gene_name, hunter_direction, hunter_log2FC = round(hunter_log2FC, 2),
              our_sig_lane, our_norm_log2FC = round(our_norm_log2FC, 2), verdict)])

# =============================================================================
# STEP 4: Our genes -> Hunter's table
# =============================================================================
# "Not flagged" and "not expressed in their LCLs" are different answers, so
# presence in their table is reported separately from significance.

hunter_chr21 <- unique(hunter[chr == "chr21", .(Gene_name, h_baseMean = baseMean,
                                                h_log2FC = log2FoldChange, h_padj = padj)],
                       by = "Gene_name")
rev <- merge(ours[, .(Gene_name, our_direction, our_sig_lane = sig_lane, our_tier = tier,
                      our_norm_log2FC = norm_log2FC, our_eqtl_lane = eqtl_lane)],
             hunter_chr21, by = "Gene_name", all.x = TRUE)
rev <- merge(rev, sig[, .(Gene_name, hunter_direction)], by = "Gene_name", all.x = TRUE)
rev[, verdict := fcase(
  is.na(h_baseMean),                            "absent from Hunter's table",
  is.na(hunter_direction),                      "in their table, not deviating for them",
  hunter_direction == our_direction,            "deviating in both, same direction",
  default =                                     "deviating in both, OPPOSITE direction")]
setorder(rev, our_direction, our_norm_log2FC)
fwrite(rev, run$table("hunter_comparison_reverse.csv"))

cat("\nStep 4: whether our deviating genes were flagged by Hunter\n")
print(rev[, .(Gene_name, our_direction, our_norm_log2FC = round(our_norm_log2FC, 2),
              h_log2FC = round(h_log2FC, 2), h_padj = signif(h_padj, 2), verdict)])

# =============================================================================
# STEP 5: Verdict
# =============================================================================

shared <- cmp[!is.na(our_direction), Gene_name]
same   <- cmp[verdict == "deviating in both, same direction", Gene_name]
opp    <- cmp[verdict == "deviating in both, OPPOSITE direction", Gene_name]

cat(sprintf("\n=== Summary [run: %s] ===\n", run$name))
cat(sprintf("  Hunter lower-than-expected (their result): %s\n",
            paste(sig[hunter_direction == "lower", Gene_name], collapse = ", ")))
cat(sprintf("  Hunter higher-than-expected (not discussed in the paper): %s\n",
            paste(sig[hunter_direction == "higher", Gene_name], collapse = ", ")))
cat(sprintf("  Deviating in both analyses: %d%s\n", length(shared),
            if (length(shared)) paste0(" (", paste(shared, collapse = ", "), ")") else ""))
cat(sprintf("    same direction:     %d%s\n", length(same),
            if (length(same)) paste0(" (", paste(same, collapse = ", "), ")") else ""))
cat(sprintf("    opposite direction: %d%s\n", length(opp),
            if (length(opp)) paste0(" (", paste(opp, collapse = ", "), ")") else ""))
cat(sprintf("  Hunter genes absent from our chr21 target set: %s\n",
            paste(cmp[verdict == "not in our chr21 target set", Gene_name], collapse = ", ")))
cat(sprintf("  Our deviating genes absent from Hunter's table: %s\n",
            paste(rev[verdict == "absent from Hunter's table", Gene_name], collapse = ", ")))
cat("\n  Wrote", run$table("hunter_comparison.csv"), "\n")
cat("  Wrote", run$table("hunter_comparison_reverse.csv"), "\n")
cat("\n=== Comparison complete ===\n")
