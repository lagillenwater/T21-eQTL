# 17_distal_inflation.R
#
# Purpose: Does the T21 cohort carry genotype-expression structure (ancestry,
#          relatedness, sequencing or sample batch) that would inflate the
#          within-T21 eQTL tests? The decoy control of script 03 is too small
#          to measure calibration (50 tests on 10 genes' expression), so this
#          measures it directly: the run's covariate-adjusted T21 expression
#          of every assessable chr21 gene (not repeat-flagged, not low
#          expression) is regressed on every common chr21 variant at least
#          decoy_min_distance (5 Mb) from the gene's TSS, one ordinary
#          regression per pair with no permutation, as in fit_variants(). At
#          that distance a variant cannot tag the gene's own cis locus, so
#          the p-values should be uniform; the genomic inflation factor
#          lambda and the QQ curve measure any departure.
#
#   Two references read the result:
#   - chr22 genes (target biotypes, baseMean >= low_expr_basemean, not
#     blacklisted) against the same variants: every pair is on another
#     chromosome. Ancestry or relatedness inflates chr22 genes as much as
#     chr21 genes; a chr21-only excess points at something chr21-specific
#     (residual mosaicism or chr21 copy-number signal in the genotype calls).
#   - A subject permutation applied jointly to all genes: it keeps
#     co-expression and LD, so its spread is the sampling spread of lambda
#     under these correlated tests. Observed lambda outside it is structure.
#   Expression is also rank inverse-normal transformed as a sensitivity: an
#   excess that disappears under the transform comes from outlying
#   expression, not from structure. Short-arm and pericentromeric variants
#   (< 13 Mb) are reported separately as well: GRCh38 falsely duplicates
#   parts of 21q22.3 there, so their coordinates can put a q-arm gene's own
#   cis signal more than 5 Mb away.
#
# Inputs (per run):
#   - processed/analysis_cohort.csv, expression_adjusted.csv
#   - tables/chr21_lane_assignments.csv, deseq2_all_genes_ploidy_normalized.csv
#   - data/chr21_ds_PASS.csv (streamed once to the cache below)
#   - data/gencode.v44.basic.annotation.gtf.gz (gene TSS)
#   - data/processed/blacklisted_genes.csv
# Outputs:
#   - data/processed/t21_chr21_dosage_maf02.csv   cache, run-independent:
#     ALT dosage of every biallelic PASS variant at file-wide MAF >= 0.02
#   - tables/distal_inflation_summary.csv  lambda per gene set and MAF bin,
#     with the permutation reference and the inverse-normal sensitivity
#   - tables/distal_inflation_genes.csv    per gene: lambda, its own
#     permutation range, smallest distal p
#   - tables/distal_inflation_gene_excess.csv  genes above their own
#     permutation 95th percentile, against the permutation count
#   - tables/distal_inflation_covariates.csv   runs with covariates: each
#     covariate regressed on every common variant (does it track genotype?)
#   - tables/distal_inflation_qq.csv       QQ curves with permutation envelope
#   - tables/distal_inflation_top_pairs.csv distal pairs at p < 1e-6
#   - figures/distal_inflation.{pdf,png}
#
# Usage: Rscript scripts/17_distal_inflation.R --run adjusted
#   T21_N_PERM_LAMBDA (default 200) subject permutations;
#   T21_CORES (default 6) parallel workers for them.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(parallel)
})
source("scripts/lib/cohort.R")
source("scripts/lib/biotypes.R")
source("scripts/lib/inflation.R")
source("scripts/lib/run.R"); run <- load_run()

N_PERM   <- as.integer(Sys.getenv("T21_N_PERM_LAMBDA", "200"))
N_CORES  <- as.integer(Sys.getenv("T21_CORES", "6"))
SEED     <- 1717L
th       <- run$thresholds
MIN_DIST <- th$decoy_min_distance   # same rule as the decoy sets: no LD with the gene's locus
MIN_MAF  <- th$decoy_min_maf        # "common": the decoy MAF floor, here in the run's T21 subjects
MIN_CALL <- 0.95
MAF_BREAKS <- c(0.05, 0.10, 0.20, 0.50)
HIT_X    <- 6                       # record distal pairs at p < 1e-6
PARM_END <- 13e6                    # end of the chr21 short arm and centromere (GRCh38 cen 10.9-12.9 Mb)
CONTROL_CHR  <- "chr22"
CONTROL_SET  <- "chr22 reference"
DOSAGE_CACHE <- "data/processed/t21_chr21_dosage_maf02.csv"
CACHE_MAF    <- 0.02
GTF <- "data/gencode.v44.basic.annotation.gtf.gz"
COL_CHR <- c(chr21 = "#B2182B", chr22 = "#2166AC")

cat(sprintf("=== T21-eQTL: distal-variant inflation [run: %s] ===\n\n", run$name))

# ---- Genotype cache (run-independent) ----------------------------------------
if (!file.exists(DOSAGE_CACHE)) {
  cat("  Streaming data/chr21_ds_PASS.csv to the dosage cache (~3 min) ...\n")
  part <- paste0(DOSAGE_CACHE, ".part")
  status <- system2("awk", c("-v", sprintf("minmaf=%s", CACHE_MAF), "-f",
                             "scripts/lib/t21_dosage_stream.awk", "data/chr21_ds_PASS.csv"),
                    stdout = part)
  if (status != 0) stop("awk stream failed")
  file.rename(part, DOSAGE_CACHE)
}

# ---- Subjects: the run's T21 cohort with expression and genotypes ------------
E <- run$expression()
t21 <- run$cohort()[Karyotype == "T21" & (is.na(excluded_reason) | excluded_reason == "")]
t21[, subject_id := subject_id_from_labid(LabID)]
t21 <- t21[LabID %in% colnames(E)]

dose <- fread(DOSAGE_CACHE, showProgress = FALSE)
# A position carried by two rows is a split multi-allelic site; drop both.
dose <- dose[!POS %in% POS[duplicated(POS)]]
G <- dosage_matrix(dose, t21$subject_id)
has_geno <- rowSums(!is.na(G)) > 0
t21 <- t21[has_geno]; G <- G[has_geno, , drop = FALSE]

vs <- dosage_stats(G)
keep_v <- vs$call_rate >= MIN_CALL & vs$maf >= MIN_MAF
G   <- G[, keep_v, drop = FALSE]
vs  <- vs[keep_v]
pos <- dose$POS[keep_v]
MAF_LAB <- sprintf("%.2f-%.2f", head(MAF_BREAKS, -1), MAF_BREAKS[-1])
maf_bin <- findInterval(vs$maf, MAF_BREAKS, rightmost.closed = TRUE)
# Variant category = arm x MAF bin. Short-arm and pericentromeric positions
# (< PARM_END) include GRCh38's false duplications of 21q22.3 (Behera et al.
# 2023), where a "distal" coordinate is really a mis-placed copy of a q-arm
# gene's own cis region; the q-arm-only set is the sensitivity that removes them.
arm <- 1L + (pos >= PARM_END)
variant_bin <- (arm - 1L) * length(MAF_LAB) + maf_bin
VBINS <- c(list("all" = seq_len(2L * length(MAF_LAB)),
                "q-arm only" = length(MAF_LAB) + seq_along(MAF_LAB)),
           setNames(lapply(seq_along(MAF_LAB), function(k) c(k, k + length(MAF_LAB))), MAF_LAB))
n_var <- vapply(VBINS, function(b) sum(variant_bin %in% b), integer(1))
cat(sprintf("  T21 subjects: %d; cache variants: %d; common (MAF >= %.2f, call rate >= %.2f): %d (%d on the q-arm)\n",
            nrow(t21), nrow(dose), MIN_MAF, MIN_CALL, ncol(G), n_var[["q-arm only"]]))

# ---- Genes: assessable chr21 and the chr22 reference -------------------------
lanes <- fread(run$table("chr21_lane_assignments.csv"))
chr21 <- lanes[!sig_lane %in% c("High_repeats", "Low_expression", "Not_assessable"),
               .(Gene_name, ensembl_stable, baseMean, sig_lane)]
chr21[, gene_set := fcase(sig_lane %in% c("DE_high", "DE_low"), "Deviating",
                          sig_lane == "Expected_dosage", "Expected dosage",
                          default = "Past magnitude, padj >= alpha")]
de <- fread(run$table("deseq2_all_genes_ploidy_normalized.csv"))
blacklist <- fread("data/processed/blacklisted_genes.csv")[in_blacklist == TRUE, Gene_name]
dup_names <- rownames(E)[duplicated(rownames(E))]
chr22 <- de[Chr == CONTROL_CHR & Gene_type %in% TARGET_BIOTYPES &
              baseMean >= th$low_expr_basemean & !Gene_name %in% c(blacklist, dup_names) &
              Gene_name %in% rownames(E),
            .(Gene_name, ensembl_stable = sub("\\..*$", "", EnsemblID), baseMean,
              sig_lane = NA_character_, gene_set = CONTROL_SET)]
stopifnot(!anyDuplicated(chr21$Gene_name), all(chr21$Gene_name %in% rownames(E)),
          !any(chr21$Gene_name %in% dup_names))

genes <- rbind(chr21[, chromosome := "chr21"], chr22[, chromosome := CONTROL_CHR])
tss <- gene_tss_gtf(GTF, c("chr21", CONTROL_CHR))
genes <- merge(genes, tss[, .(ensembl_stable, tss)], by = "ensembl_stable", all.x = TRUE, sort = FALSE)
if (genes[chromosome == "chr21", anyNA(tss)]) stop("chr21 gene without a GENCODE TSS")
GROUPS <- c("Deviating", "Expected dosage", "Past magnitude, padj >= alpha", CONTROL_SET)
genes[, group := match(gene_set, GROUPS)]
setorder(genes, group, Gene_name)
gene_tss <- genes[, fifelse(chromosome == "chr21", as.numeric(tss), NA_real_)]
cat(sprintf("  Genes: %s\n", paste(sprintf("%s %d", GROUPS, tabulate(genes$group, 4)), collapse = "; ")))

# Gene sets reported: all assessable chr21 genes, each chr21 group, chr22. A
# set with no genes in this run (the baseline has no past-magnitude genes) is
# dropped rather than reported as NA.
SETS <- list("chr21 assessable" = 1:3, "chr21 deviating" = 1L, "chr21 expected dosage" = 2L,
             "chr21 past magnitude, padj >= alpha" = 3L, "chr22 reference" = 4L)
SETS <- SETS[vapply(SETS, function(g) any(genes$group %in% g), logical(1))]
LAMBDA_GRID <- CJ(gene_set = names(SETS), variant_set = names(VBINS), sorted = FALSE)
QQ_GRID <- CJ(gene_set = c("chr21 assessable", CONTROL_SET), variant_set = c("all", "q-arm only"),
              sorted = FALSE)
set_lambdas <- function(counts) LAMBDA_GRID[, {
  sub <- counts[, VBINS[[variant_set]], SETS[[gene_set]], drop = FALSE]
  .(n_tests = sum(sub), lambda = hist_lambda(sub, bins))
}, by = .(gene_set, variant_set)]
set_qq <- function(counts) QQ_GRID[, hist_qq(counts[, VBINS[[variant_set]], SETS[[gene_set]], drop = FALSE],
                                             bins), by = .(gene_set, variant_set)]

Y  <- t(E[genes$Gene_name, t21$LabID, drop = FALSE])
stopifnot(!anyNA(Y))
n  <- nrow(Y); df <- n - 2L
Ys <- standardise_cols(Y)
Gs <- standardise_cols(G)
bins <- log10p_bins(df)

scan_args <- list(Gs = Gs, gene_tss = gene_tss, pos = pos, gene_group = genes$group,
                  variant_bin = variant_bin, min_dist = MIN_DIST, bins = bins,
                  n_groups = length(GROUPS), n_maf = 2L * length(MAF_LAB))

# ---- Observed scan, and the inverse-normal sensitivity -----------------------
t0 <- Sys.time()
obs <- do.call(distal_scan, c(list(Ys = Ys, hit_x = HIT_X), scan_args))
cat(sprintf("  Observed scan: %s distal tests (%.1f min)\n",
            format(sum(obs$counts), big.mark = ","), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
int <- do.call(distal_scan, c(list(Ys = standardise_cols(inverse_normal_cols(Y))), scan_args))

# ---- Permutation reference ---------------------------------------------------
cat(sprintf("  Permutation reference: %d subject permutations on %d workers ...\n", N_PERM, N_CORES))
t0 <- Sys.time()
perm <- mclapply(seq_len(N_PERM), function(i) {
  setDTthreads(1L)
  set.seed(SEED + i)
  s <- do.call(distal_scan, c(list(Ys = Ys[sample.int(n), , drop = FALSE]), scan_args))
  list(lambda = set_lambdas(s$counts)[, perm := i], qq = set_qq(s$counts)[, perm := i],
       gene_median_r = s$gene_median_r)
}, mc.cores = N_CORES, mc.preschedule = FALSE)
failed <- vapply(perm, inherits, logical(1), "try-error")
if (any(failed)) stop(sprintf("%d permutation workers failed: %s", sum(failed), perm[[which(failed)[1]]]))
cat(sprintf("  done (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

perm_lambda <- rbindlist(lapply(perm, `[[`, "lambda"))
perm_qq     <- rbindlist(lapply(perm, `[[`, "qq"))
perm_gene   <- vapply(perm, `[[`, numeric(nrow(genes)), "gene_median_r")   # genes x perms

# ---- Summary table -----------------------------------------------------------
KEY <- c("gene_set", "variant_set")
summary_tbl <- set_lambdas(obs$counts)
summary_tbl <- merge(summary_tbl, set_lambdas(int$counts)[, .(gene_set, variant_set, lambda_inverse_normal = lambda)],
                     by = KEY, sort = FALSE)
perm_ref <- merge(perm_lambda, summary_tbl[, .(gene_set, variant_set, observed = lambda)], by = KEY)[
  , .(perm_median = median(lambda), perm_lo = quantile(lambda, 0.025),
      perm_hi = quantile(lambda, 0.975), perm_p = (1 + sum(lambda >= observed)) / (N_PERM + 1)),
  by = KEY]
summary_tbl <- merge(summary_tbl, perm_ref, by = KEY, sort = FALSE)
summary_tbl[, n_genes := vapply(gene_set, function(s) sum(genes$group %in% SETS[[s]]), integer(1))]
summary_tbl[, n_variants := n_var[variant_set]]
summary_tbl[, `:=`(n_subjects = n, min_distance_bp = MIN_DIST, n_perm = N_PERM)]
setcolorder(summary_tbl, c(KEY, "n_genes", "n_variants", "n_tests", "lambda",
                           "perm_median", "perm_lo", "perm_hi", "perm_p", "lambda_inverse_normal"))
num <- c("lambda", "perm_median", "perm_lo", "perm_hi", "lambda_inverse_normal")
summary_tbl[, (num) := lapply(.SD, round, 4), .SDcols = num]
summary_tbl[, perm_p := signif(perm_p, 3)]
fwrite(summary_tbl, run$table("distal_inflation_summary.csv"))

# ---- Per-gene table ----------------------------------------------------------
gene_lambda_obs  <- lambda_from_abs_r(obs$gene_median_r, df)
gene_lambda_perm <- matrix(lambda_from_abs_r(perm_gene, df), nrow(genes))
thr95 <- apply(gene_lambda_perm, 1, quantile, 0.95)
genes[, `:=`(
  n_tests = obs$gene_n,
  lambda = round(gene_lambda_obs, 4),
  lambda_inverse_normal = round(lambda_from_abs_r(int$gene_median_r, df), 4),
  perm_lo = round(apply(gene_lambda_perm, 1, quantile, 0.025), 4),
  perm_hi = round(apply(gene_lambda_perm, 1, quantile, 0.975), 4),
  perm_p  = signif((1 + rowSums(gene_lambda_perm >= gene_lambda_obs)) / (N_PERM + 1), 3),
  above_own_95 = gene_lambda_obs > thr95,
  min_distal_p = signif(p_from_abs_r(obs$gene_max_r, df), 3))]
gene_out <- genes[, .(Gene_name, ensembl_stable, chromosome, gene_set, sig_lane, baseMean = round(baseMean, 1),
                      tss, n_tests, lambda, perm_lo, perm_hi, perm_p, above_own_95,
                      lambda_inverse_normal, min_distal_p)]
setorder(gene_out, -lambda)
fwrite(gene_out, run$table("distal_inflation_genes.csv"))

# ---- Excess of genes above their own permutation range -----------------------
# Each gene is compared with its own 95th permutation percentile; the null
# spread of the count comes from applying the same thresholds to every
# permutation, which keeps co-expression (a binomial would not).
above_perm <- gene_lambda_perm > thr95
excess <- rbindlist(lapply(names(SETS), function(s) {
  w <- genes$group %in% SETS[[s]]
  cnt <- colSums(above_perm[w, , drop = FALSE])
  k <- sum(genes$above_own_95[w])
  data.table(gene_set = s, n_genes = sum(w), n_above_own_95 = k,
             perm_median = median(cnt), perm_hi = unname(quantile(cnt, 0.95)),
             perm_p = signif((1 + sum(cnt >= k)) / (N_PERM + 1), 3))
}))
fwrite(excess, run$table("distal_inflation_gene_excess.csv"))

# ---- Do the run's covariates themselves track chr21 genotype? ----------------
# Each covariate is regressed on every common variant (no locus, so no distance
# rule), with the same subject-permutation reference. A covariate associated
# with genotype across the whole chromosome shares the structure (ancestry,
# recruitment) that would inflate the eQTL tests; adjusting for it removes that
# structure by proxy. Runs without covariates skip this.
X <- if (length(run$covariates)) {
  model.matrix(~ ., t21[, run$covariates, with = FALSE])[, -1, drop = FALSE]
} else NULL
if (!is.null(run$composition)) {
  comp <- fread(run$processed("composition_fractions.csv"))
  X <- cbind(X, as.matrix(comp[match(t21$LabID, comp$LabID), !"LabID"]))
}
cov_out <- NULL
if (!is.null(X)) {
  stopifnot(nrow(X) == n, !anyNA(X))
  Xs <- standardise_cols(X)
  cov_lambda <- function(M) lambda_from_abs_r(
    distal_scan(M, Gs, rep(NA_real_, ncol(M)), pos, rep(1L, ncol(M)), variant_bin, 0, bins,
                n_groups = 1L, n_maf = 2L * length(MAF_LAB))$gene_median_r, df)
  cov_obs  <- cov_lambda(Xs)
  cov_perm <- simplify2array(mclapply(seq_len(N_PERM), function(i) {
    set.seed(SEED + 10000L + i)
    cov_lambda(Xs[sample.int(n), , drop = FALSE])
  }, mc.cores = N_CORES))                                     # covariates x perms
  cov_thr95 <- apply(cov_perm, 1, quantile, 0.95)
  cov_out <- data.table(covariate = colnames(X), lambda = round(cov_obs, 4),
                        perm_lo = round(apply(cov_perm, 1, quantile, 0.025), 4),
                        perm_hi = round(apply(cov_perm, 1, quantile, 0.975), 4),
                        perm_p = signif((1 + rowSums(cov_perm >= cov_obs)) / (N_PERM + 1), 3),
                        above_own_95 = cov_obs > cov_thr95)
  setorder(cov_out, -lambda)
  fwrite(cov_out, run$table("distal_inflation_covariates.csv"))
  cnt <- colSums(cov_perm > cov_thr95)
  excess <- rbind(excess, data.table(
    gene_set = "run covariates (not genes)", n_genes = ncol(X), n_above_own_95 = sum(cov_out$above_own_95),
    perm_median = median(cnt), perm_hi = unname(quantile(cnt, 0.95)),
    perm_p = signif((1 + sum(cnt >= sum(cov_out$above_own_95))) / (N_PERM + 1), 3)))
  fwrite(excess, run$table("distal_inflation_gene_excess.csv"))
}

# ---- QQ table ----------------------------------------------------------------
qq_obs <- set_qq(obs$counts)
qq_env <- perm_qq[, .(perm_lo = quantile(observed, 0.025), perm_hi = quantile(observed, 0.975)),
                  by = c(KEY, "rank")]
qq_tbl <- merge(qq_obs, qq_env, by = c(KEY, "rank"), sort = FALSE)
fwrite(qq_tbl, run$table("distal_inflation_qq.csv"))

# ---- Top distal pairs --------------------------------------------------------
hits <- if (is.null(obs$hits)) data.table(gene = integer(0), variant = integer(0), r = numeric(0)) else obs$hits
top <- hits[, .(Gene_name = genes$Gene_name[gene], gene_set = genes$gene_set[gene],
                variant_id = colnames(G)[variant], q_arm = pos[variant] >= PARM_END,
                maf = round(vs$maf[variant], 3),
                distance_mb = round(abs(pos[variant] - gene_tss[gene]) / 1e6, 2),   # NA off chr21
                r = round(r, 4), p = signif(p_from_abs_r(r, df), 3))]
setorder(top, p)
fwrite(top, run$table("distal_inflation_top_pairs.csv"))

# ---- Figure ------------------------------------------------------------------
qq_tbl <- merge(qq_tbl, summary_tbl[, .(gene_set, variant_set, n_genes, lambda, perm_lo_l = perm_lo,
                                        perm_hi_l = perm_hi)], by = KEY)
qq_tbl[, curve := sprintf("%s (%d genes), %s variants: lambda %.3f (permutation 95%% %.3f-%.3f)",
                          gene_set, n_genes, sub(" only", "", variant_set), lambda, perm_lo_l, perm_hi_l)]
curves <- unique(qq_tbl[, .(curve, gene_set, variant_set)])
col_curve <- setNames(COL_CHR[fifelse(curves$gene_set == CONTROL_SET, "chr22", "chr21")], curves$curve)
lty_curve <- setNames(fifelse(curves$variant_set == "all", "solid", "22"), curves$curve)
# The axis is set by the q-arm curves and the envelopes; the short-arm
# false-duplication tail of the all-variant chr21 curve runs off the top.
y_cap <- ceiling(max(qq_tbl[variant_set == "q-arm only" | gene_set == CONTROL_SET, observed],
                     qq_tbl$perm_hi, qq_tbl$expected)) + 1
off_top <- qq_tbl[observed > y_cap, max(observed)]
pA <- ggplot(qq_tbl, aes(expected, observed, colour = curve, linetype = curve)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_ribbon(data = qq_tbl[variant_set == "all"], aes(ymin = perm_lo, ymax = perm_hi, fill = curve),
              colour = NA, alpha = 0.18, show.legend = FALSE) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = col_curve, name = NULL) +
  scale_fill_manual(values = col_curve, guide = "none") +
  scale_linetype_manual(values = lty_curve, name = NULL) +
  coord_cartesian(xlim = c(0, y_cap), ylim = c(0, y_cap), expand = FALSE) +
  labs(x = expression(Expected ~ -log[10] ~ italic(p)), y = expression(Observed ~ -log[10] ~ italic(p)),
       title = "A  Distal genotype-expression tests, T21 subjects",
       subtitle = sprintf("Common chr21 variants >= %.0f Mb from each TSS; shaded: permutation 95%% envelope",
                          MIN_DIST / 1e6)) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.key.width = unit(1.4, "lines"),
        panel.grid.minor = element_blank(), aspect.ratio = 1)
if (length(off_top) && is.finite(off_top)) {
  pA <- pA + annotate("text", x = 0.3, y = y_cap - 0.4, hjust = 0, vjust = 1, size = 2.5, colour = "grey25",
                      label = sprintf("chr21, all variants: continues to %.0f\n(short-arm false-duplication variants)",
                                      off_top))
}

null_band <- quantile(gene_lambda_perm, c(0.025, 0.975))
gene_out[, set_f := factor(gene_set, levels = GROUPS[GROUPS %in% gene_set])]
# Label the deviating genes and the eight highest lambdas elsewhere.
gene_out[, label := fifelse(gene_set == "Deviating" | rank(-lambda, ties.method = "first") <= 8,
                            Gene_name, "")]
pB <- ggplot(gene_out, aes(set_f, lambda, colour = chromosome)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = null_band[1], ymax = null_band[2],
           fill = "grey80", alpha = 0.5) +
  geom_hline(yintercept = 1, colour = "grey55", linewidth = 0.3) +
  geom_point(position = position_jitter(width = 0.25, height = 0, seed = SEED), size = 1.2, alpha = 0.8) +
  geom_text_repel(aes(label = label), size = 2.3, colour = "grey15", max.overlaps = 30,
                  position = position_jitter(width = 0.25, height = 0, seed = SEED),
                  segment.size = 0.2, min.segment.length = 0) +
  scale_colour_manual(values = COL_CHR, guide = "none") +
  scale_x_discrete(labels = function(x) sub(", ", ",\n", x)) +
  labs(x = NULL, y = "Per-gene lambda",
       title = "B  Per-gene lambda",
       subtitle = "Grey band: 95% of per-gene lambda under subject permutation") +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())

var_tbl <- summary_tbl[gene_set %in% QQ_GRID$gene_set]
var_tbl[, var_f := factor(variant_set, levels = names(VBINS))]
var_tbl[, chromosome := fifelse(gene_set == CONTROL_SET, "chr22", "chr21")]
dodge <- position_dodge(width = 0.5)
pC <- ggplot(var_tbl, aes(var_f, lambda, colour = chromosome, group = chromosome)) +
  geom_hline(yintercept = 1, colour = "grey55", linewidth = 0.3) +
  geom_linerange(aes(ymin = perm_lo, ymax = perm_hi), position = dodge, linewidth = 2.2,
                 colour = "grey80") +
  geom_point(position = dodge, size = 2) +
  geom_point(aes(y = lambda_inverse_normal), position = dodge, size = 2, shape = 1) +
  scale_colour_manual(values = COL_CHR, labels = c(chr21 = "chr21 assessable", chr22 = CONTROL_SET),
                      name = NULL) +
  labs(x = "Variants (T21 minor allele frequency bins use both arms)", y = "lambda",
       title = "C  lambda by variant set",
       subtitle = "Filled: adjusted expression; open: rank inverse-normal; grey: permutation 95%") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank())

fig <- pA + (pB / pC) + plot_layout(widths = c(1.1, 1))
fig_h <- 7.5
if (!is.null(cov_out)) {
  cov_out[, covariate_f := factor(covariate, levels = rev(covariate))]
  pD <- ggplot(cov_out, aes(lambda, covariate_f)) +
    geom_vline(xintercept = 1, colour = "grey55", linewidth = 0.3) +
    geom_linerange(aes(xmin = perm_lo, xmax = perm_hi), linewidth = 2.2, colour = "grey80") +
    geom_point(size = 1.8, colour = "grey15") +
    labs(x = "lambda (covariate ~ each common chr21 variant)", y = NULL,
         title = "D  Run covariates against chr21 genotype",
         subtitle = "Grey: permutation 95%") +
    theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
  fig <- (pA | pB) / (pC | pD)
  fig_h <- 12
}
stem <- run$figure("distal_inflation")
ggsave(paste0(stem, ".pdf"), fig, width = 11, height = fig_h, device = cairo_pdf)   # cairo: Greek in cell-type names
ggsave(paste0(stem, ".png"), fig, width = 11, height = fig_h, dpi = 200)

cat("\n"); print(summary_tbl[variant_set %in% c("all", "q-arm only"),
                             .(gene_set, variant_set, n_genes, n_tests, lambda, perm_lo, perm_hi,
                               perm_p, lambda_inverse_normal)])
cat("\n"); print(excess)
cat(sprintf("\n  Distal pairs at p < 1e-%d: %d (%d genes; %d with a short-arm variant)\n", HIT_X,
            nrow(top), uniqueN(top$Gene_name), sum(!top$q_arm)))
writeLines(capture.output(sessionInfo()), run$table("distal_inflation_session_info.txt"))
cat("\nDone.\n")
