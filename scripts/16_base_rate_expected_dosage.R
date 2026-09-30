# 16_base_rate_expected_dosage.R
#
# Purpose: Base rate for the cis-eQTL detection result (manuscript review,
#          concern 2). "7 of 10 deviating genes have a detected cis-eQTL"
#          can only be read against how often the identical test detects one
#          for a chr21 gene that sits at the expected dosage. This script runs
#          the same gene-level permutation test (same variant universe, same
#          subjects, same expression artifact, same thresholds) on every
#          Expected-dosage gene of the run that has GTEx cis variants at the
#          pval cut, and compares the detection rates.
#
#   Results are appended per chunk of genes, so a partial run still leaves a
#   usable table; BH is applied over whatever finished.
#
# Inputs (per run):
#   - tables/chr21_lane_assignments.csv, eqtl_gene_level_perm.csv
#   - processed/analysis_cohort.csv, expression_adjusted.csv
#   - data/GTEx allpairs chr21 parquet, data/chr21_ds_PASS.csv
# Outputs (per run):
#   - tables/base_rate_expected_dosage.csv        one row per tested gene
#   - tables/base_rate_summary.csv                 rates: deviating vs expected dosage
#   - processed/base_rate_genotypes.csv            T21 dosage at the pulled variants
#
# Usage: Rscript scripts/16_base_rate_expected_dosage.R --run adjusted
#   T21_N_PERM (default 1000) sets the permutations per gene.

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})
source("scripts/lib/cohort.R")
source("scripts/lib/eqtl_fit.R")
source("scripts/lib/eqtl_controls.R")
source("scripts/lib/run.R"); run <- load_run()

N_PERM   <- as.integer(Sys.getenv("T21_N_PERM", "1000"))
CHUNK    <- 10L
th       <- run$thresholds
ALLPAIRS <- paste0("data/GTEx_Analysis_v10_QTLs_GTEx_Analysis_v10_eQTL",
                   "_all_associations_Whole_Blood.v10.allpairs.chr21.parquet")

cat(sprintf("=== T21-eQTL: base rate on Expected-dosage genes [run: %s] ===\n\n", run$name))

# ---- Genes and their GTEx cis variants at the cut ---------------------------
lanes <- fread(run$table("chr21_lane_assignments.csv"))
lanes[, ensembl_stable := sub("\\..*$", "", EnsemblID)]
expected <- lanes[sig_lane == "Expected_dosage", .(Gene_name, ensembl_stable, baseMean, norm_log2FC)]
gtex <- as.data.table(read_parquet(ALLPAIRS))
gtex[, ensembl_stable := sub("\\..*$", "", gene_id)]
tv <- gtex[ensembl_stable %in% expected$ensembl_stable & startsWith(variant_id, "chr21_") &
             !is.na(pval_nominal) & pval_nominal <= th$gtex_pval_keep,
           .(ensembl_stable, variant_id, pval_nominal)]
parsed <- tstrsplit(tv$variant_id, "_", fixed = TRUE)
tv[, `:=`(POS = as.integer(parsed[[2]]), REF = parsed[[3]], ALT = parsed[[4]])]
tv <- merge(tv, expected[, .(ensembl_stable, Gene_name)], by = "ensembl_stable")
cat(sprintf("  Expected-dosage genes: %d; with GTEx variants at p <= %.0e: %d; variant rows: %d\n",
            nrow(expected), th$gtex_pval_keep, uniqueN(tv$Gene_name), nrow(tv)))

# ---- T21 genotypes at those positions (awk stream, T21 file only) -----------
cohort <- run$cohort()
t21 <- cohort[Karyotype == "T21" & (is.na(excluded_reason) | excluded_reason == "")]
t21[, subject_id := subject_id_from_labid(LabID)]
geno_out <- run$processed("base_rate_genotypes.csv")
if (file.exists(geno_out)) {
  cat("  Genotype extract exists; reusing.\n")
  geno <- fread(geno_out)
} else {
  pos_tmp <- tempfile(fileext = ".txt"); writeLines(as.character(unique(tv$POS)), pos_tmp)
  awk_cmd <- sprintf(paste0("awk -F',' 'BEGIN{while((getline l < \"%s\")>0) p[l]=1} ",
                            "NR==1 || ($2 in p)' %s"), pos_tmp, "data/chr21_ds_PASS.csv")
  cat("  Streaming data/chr21_ds_PASS.csv ...\n")
  dt <- fread(cmd = awk_cmd, sep = ",", header = TRUE, showProgress = FALSE)
  meta_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  geno_cols <- setdiff(names(dt), meta_cols)
  keep <- geno_cols[sub("[A-Z][0-9]*$", "", geno_cols) %in% t21$subject_id]
  dt <- dt[!grepl(",", ALT), c("POS", "REF", "ALT", keep), with = FALSE]
  dt <- merge(dt, unique(tv[, .(variant_id, POS, REF, ALT)]), by = c("POS", "REF", "ALT"))
  long <- melt(dt, id.vars = c("variant_id", "POS", "REF", "ALT"), measure.vars = keep,
               variable.name = "lab_id", value.name = "geno_field", variable.factor = FALSE)
  gt <- sub(":.*$", "", long$geno_field)
  long[, alt_dosage := fifelse(grepl("\\.", gt), NA_real_, nchar(gsub("[^1-9]", "", gt)))]
  long[, subject_id := sub("[A-Z][0-9]*$", "", lab_id)]
  geno <- long[!is.na(alt_dosage), .(alt_dosage = mean(alt_dosage)), by = .(variant_id, subject_id)]
  fwrite(geno, geno_out)
}
cat(sprintf("  Genotyped variants: %d; T21 subjects: %d\n", uniqueN(geno$variant_id), uniqueN(geno$subject_id)))

gwide <- dcast(geno, subject_id ~ variant_id, value.var = "alt_dosage")
subj  <- gwide$subject_id
G_all <- as.matrix(gwide[, -1])
E     <- run$expression()
expr_of <- function(g) expr_from_matrix(g, subj, E, t21)
variants_of <- split(tv$variant_id, tv$Gene_name)
genes <- names(variants_of)
genes <- genes[genes %in% rownames(E)]

# ---- The identical test, chunked so a partial run leaves results ------------
out_path <- run$table("base_rate_expected_dosage.csv")
if (file.exists(out_path)) file.remove(out_path)
chunks <- split(genes, ceiling(seq_along(genes) / CHUNK))
t0 <- Sys.time()
for (i in seq_along(chunks)) {
  res <- gene_level_tests(chunks[[i]], variants_of, G_all, expr_of,
                          n_perm = N_PERM, seed_base = 5026L + (i - 1L) * CHUNK, fdr = th$fdr_gene)
  res[, c("q_gene_bh", "detected") := NULL]
  fwrite(res, out_path, append = file.exists(out_path))
  cat(sprintf("  chunk %d/%d done (%d genes, %.1f min elapsed)\n", i, length(chunks),
              length(chunks[[i]]), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# ---- Rates ------------------------------------------------------------------
res <- fread(out_path)
res <- merge(res, expected[, .(Gene_name, baseMean, norm_log2FC)], by = "Gene_name")
res[, q_gene_bh := p.adjust(p_gene_perm, "BH")]
res[, detected := !is.na(q_gene_bh) & q_gene_bh < th$fdr_gene]
setorder(res, p_gene_perm, na.last = TRUE)
fwrite(res, out_path)

dev <- fread(run$table("eqtl_gene_level_perm.csv"))[!is.na(p_gene_perm)]
dev <- merge(dev, lanes[, .(Gene_name, baseMean)], by = "Gene_name")
bm_range <- range(dev$baseMean)
rate <- function(d, label) data.table(
  set = label, n_tested = sum(!is.na(d$p_gene_perm)), n_detected = sum(d$detected),
  pct_detected = round(100 * mean(d$detected[!is.na(d$p_gene_perm)]), 1),
  n_nominal_p05 = sum(d$p_gene_perm < 0.05, na.rm = TRUE), median_baseMean = round(median(d$baseMean)))
summary_tbl <- rbind(
  rate(dev[, .(p_gene_perm, detected = cis_eqtl_detected, baseMean)], "deviating genes"),
  rate(res, "expected-dosage genes"),
  rate(res[baseMean >= bm_range[1] & baseMean <= bm_range[2]],
       sprintf("expected-dosage genes, baseMean %.0f-%.0f (deviating range)", bm_range[1], bm_range[2])))
# Fisher test: detection among deviating vs expected-dosage genes
ft <- fisher.test(matrix(c(summary_tbl$n_detected[1], summary_tbl$n_tested[1] - summary_tbl$n_detected[1],
                           summary_tbl$n_detected[2], summary_tbl$n_tested[2] - summary_tbl$n_detected[2]), 2))
summary_tbl[, fisher_p_vs_deviating := c(NA, ft$p.value, NA)]
fwrite(summary_tbl, run$table("base_rate_summary.csv"))
cat("\n"); print(summary_tbl)
cat("\nDone.\n")
