# S1_covariate_evidence.R
#
# Purpose: Evidence for the covariates of a run, as a supplement. Four parts,
#          each a table (S1_*.csv under the run) and a panel of one figure:
#
#   A. Covariates vs karyotype on the run's pre-exclusion roster: distribution
#      by group and a test (Wilcoxon for continuous, chi-square for
#      categorical). Shows which candidate covariates are confounded with
#      karyotype in this cohort.
#   B. Composition structure and reach (composition runs only): per cell
#      fraction, mean by karyotype, log2 ratio, Wilcoxon q, and the number of
#      expressed genes associated with that fraction at 5% FDR within
#      controls and within T21 (Spearman on unadjusted log2-CPM). A PCA of
#      the fractions is reported for description only.
#   C. Mosaic karyotype vs expression: the per-sample chr21 dosage index by
#      subtype, and a genome-wide count of genes differing between mosaic and
#      full-trisomy subjects at 5% FDR, chr21 and non-chr21 separately.
#   D. Attribution of the adjustment (needs the baseline run and a
#      composition run): for each deviating gene, the change in fold change
#      between runs and its per-covariate decomposition on the artifact scale
#      (scripts/lib/covariates.R::attribution), with bootstrap intervals.
#
#          Everything here is descriptive: association, not cause.
#
# Usage:  Rscript scripts/S1_covariate_evidence.R --run <name>
# Inputs: the run's cohort_roster.csv, composition_fractions.csv (if any),
#         expression_adjusted.csv; data/processed/count_matrix.csv and
#         sample_metadata.csv; results/runs/baseline/tables/chr21_lane_assignments.csv
#         and the run's lane table for panel D.
# Outputs: results/runs/<name>/tables/S1_*.csv,
#          results/runs/<name>/figures/S1_covariate_evidence.{pdf,png}

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("scripts/lib/run.R"); run <- load_run()
source("scripts/lib/covariates.R"); source("scripts/lib/cohort.R")
set.seed(1)
cat(sprintf("=== S1 covariate evidence [run: %s] ===\n\n", run$name))

roster <- fread(run$processed("cohort_roster.csv"))
# fwrite writes NA as an empty field; read it back as NA.
roster[excluded_reason == "", excluded_reason := NA_character_]
cohort <- roster[is.na(excluded_reason)]
count_mat <- fread("data/processed/count_matrix.csv")
labs_all <- intersect(setdiff(names(count_mat), c("EnsemblID", "Gene_name", "Chr")), roster$LabID)
M <- as.matrix(count_mat[, ..labs_all]); rownames(M) <- count_mat$Gene_name
is_chr21 <- count_mat$Chr == "chr21"
# Unadjusted log2-CPM (chr21-excluded library size) for the association counts.
L0 <- log2(t(t(M) / colSums(M[!is_chr21, ]) * 1e6) + 1)
pool <- run$partner_pool(count_mat, cohort$LabID)
k_all <- roster$Karyotype[match(labs_all, roster$LabID)]
st_all <- roster$karyotype_subtype[match(labs_all, roster$LabID)]

# ---- A: covariates vs karyotype ---------------------------------------------
cat("A: covariates vs karyotype (pre-exclusion roster)\n")
test_cov <- function(v, k) {
  ok <- !is.na(v)
  if (is.numeric(v)) {
    list(test = "wilcoxon", p = wilcox.test(v[ok] ~ k[ok])$p.value,
         t21 = sprintf("median %.1f", median(v[ok & k == "T21"])),
         control = sprintf("median %.1f", median(v[ok & k == "Control"])),
         n_missing = sum(!ok))
  } else {
    tb <- table(k[ok], v[ok])
    list(test = "chi-square", p = suppressWarnings(chisq.test(tb)$p.value),
         t21 = paste(sprintf("%s=%d", colnames(tb), tb["T21", ]), collapse = "; "),
         control = paste(sprintf("%s=%d", colnames(tb), tb["Control", ]), collapse = "; "),
         n_missing = sum(!ok))
  }
}
A <- rbindlist(lapply(c("Age_at_visit", "Sex", "BMI", "Sample_source", "Event_name"), function(cv)
  c(list(covariate = cv), test_cov(roster[[cv]], roster$Karyotype))))
A[, in_run_model := covariate %in% run$covariates]
fwrite(A, run$table("S1_covariates_vs_karyotype.csv")); print(as.data.frame(A))

# ---- B: composition structure and reach -------------------------------------
B <- NULL; pcs <- NULL
if (!is.null(run$composition)) {
  cat("\nB: composition structure and reach\n")
  Fr <- fread(run$processed("composition_fractions.csv"))
  F <- as.matrix(Fr[, -1]); rownames(F) <- Fr$LabID
  labs <- intersect(rownames(F), colnames(L0)); F <- F[labs, ]
  k <- cohort$Karyotype[match(labs, cohort$LabID)]
  assoc_count <- function(grp) {
    Lg <- L0[pool, labs[k == grp]]; Fg <- F[k == grp, ]
    r <- cor(t(Lg), Fg, method = "spearman"); n <- ncol(Lg)
    tt <- r * sqrt((n - 2) / pmax(1 - r^2, 1e-12)); p <- 2 * pt(-abs(tt), n - 2)
    colSums(apply(p, 2, p.adjust, method = "BH") < 0.05)
  }
  B <- data.table(cluster = colnames(F),
                  mean_t21 = colMeans(F[k == "T21", ]), mean_control = colMeans(F[k == "Control", ]))
  B[, log2_ratio := log2(mean_t21 / mean_control)]
  B[, q_wilcox := p.adjust(vapply(colnames(F), function(cn) wilcox.test(F[k == "T21", cn], F[k == "Control", cn])$p.value, 0), "BH")]
  B[, genes_assoc_control := assoc_count("Control")]
  B[, genes_assoc_t21 := assoc_count("T21")]
  pc <- prcomp(scale(F))
  pcs <- data.table(PC = seq_len(ncol(F)), var_expl = round(pc$sdev^2 / sum(pc$sdev^2), 3),
                    p_karyotype = apply(pc$x, 2, function(s) wilcox.test(s ~ k)$p.value),
                    top_loadings = apply(pc$rotation, 2, function(l) paste(names(sort(abs(l), decreasing = TRUE))[1:3], collapse = "; ")))
  fwrite(B, run$table("S1_composition_reach.csv")); fwrite(pcs, run$table("S1_composition_pcs.csv"))
  print(as.data.frame(B[order(-genes_assoc_t21)]))
}

# ---- C: mosaic karyotype vs expression --------------------------------------
cat("\nC: mosaic karyotype vs expression\n")
keep21 <- is_chr21 & rowMeans(M) >= 100
ref <- apply(L0[keep21, labs_all[k_all == "Control"], drop = FALSE], 1, median)
idx <- 2^apply(L0[keep21, labs_all, drop = FALSE] - ref, 2, median)
comp <- colSums(2^L0[keep21, labs_all, drop = FALSE] - 1) / mean(colSums(2^L0[keep21, labs_all[k_all == "Control"], drop = FALSE] - 1))
C <- data.table(LabID = labs_all, Karyotype = k_all, karyotype_subtype = st_all,
                chr21_index = idx, chr21_composite = comp,
                in_run_cohort = labs_all %in% cohort$LabID)
fwrite(C, run$table("S1_chr21_dosage_index.csv"))
print(as.data.frame(C[, .(n = .N, median_index = round(median(chr21_index), 2),
                          median_composite = round(median(chr21_composite), 2)), by = karyotype_subtype]))
full <- labs_all[st_all %in% c("T21", "DS_T21", "translocation_T21")]
mos  <- labs_all[st_all %in% "mosaic_T21"]
Cg <- NULL
if (length(mos) >= 3) {
  rows <- pool | is_chr21
  Lm <- L0[rows, mos, drop = FALSE]; Lf <- L0[rows, full, drop = FALSE]
  d  <- rowMeans(Lm) - rowMeans(Lf)
  s2 <- apply(Lm, 1, var) / length(mos) + apply(Lf, 1, var) / length(full)
  p  <- 2 * pnorm(-abs(d / sqrt(s2))); q <- p.adjust(p, "BH")
  Cg <- data.table(Gene_name = rownames(L0)[rows], chr21 = is_chr21[rows], diff_log2 = d, q = q)
  fwrite(Cg, run$table("S1_mosaic_vs_full_genes.csv"))
  print(as.data.frame(Cg[, .(n = .N, n_q05 = sum(q < 0.05), median_diff = round(median(diff_log2), 3)), by = chr21]))
}

# ---- D: attribution of the adjustment ---------------------------------------
D <- NULL
base_lanes <- file.path("results", "runs", "baseline", "tables", "chr21_lane_assignments.csv")
if (!is.null(run$composition) && file.exists(base_lanes) && file.exists(run$table("chr21_lane_assignments.csv"))) {
  cat("\nD: attribution of the adjustment\n")
  bl <- fread(base_lanes); al <- fread(run$table("chr21_lane_assignments.csv"))
  genes <- union(bl[sig_lane %in% c("DE_high", "DE_low"), Gene_name], al[sig_lane %in% c("DE_high", "DE_low"), Gene_name])
  genes <- intersect(genes, rownames(L0))
  X <- F
  if (length(run$covariates)) {
    cov_df <- as.data.frame(cohort[match(labs, LabID), run$covariates, with = FALSE])
    X <- cbind(model.matrix(~ ., cov_df)[, -1, drop = FALSE], F)
  }
  g <- k == "T21"
  boot_idx <- replicate(200, sample.int(length(labs), replace = TRUE))
  D <- rbindlist(lapply(genes, function(gn) {
    y <- L0[match(gn, rownames(L0)), labs]
    at <- attribution(y, g, X)
    bt <- apply(boot_idx, 2, function(i) attribution(y[i], g[i], X[i, , drop = FALSE])$contribution)
    data.table(Gene_name = gn, term = at$term, contribution = at$contribution,
               lo = apply(bt, 1, quantile, 0.025), hi = apply(bt, 1, quantile, 0.975),
               lfc_baseline = bl$norm_log2FC[match(gn, bl$Gene_name)],
               lfc_adjusted = al$norm_log2FC[match(gn, al$Gene_name)])
  }))
  fwrite(D, run$table("S1_attribution.csv"))
  top <- D[, .SD[which.max(abs(contribution))], by = Gene_name][, .(Gene_name, top_term = term, contribution = round(contribution, 3))]
  print(as.data.frame(top))
}

# ---- figure -----------------------------------------------------------------
cat("\nDrawing the figure...\n")
th <- theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 10))
pA <- ggplot(A, aes(x = reorder(covariate, -log10(p)), y = -log10(p), fill = in_run_model)) +
  geom_col() + coord_flip() + scale_fill_manual(values = c(`TRUE` = "#762A83", `FALSE` = "#BABABA"), name = "In run model") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  labs(title = "A  Covariates vs karyotype", y = "-log10 p", x = NULL) + th + theme(legend.position = "bottom")
pC <- ggplot(C[Karyotype == "T21"], aes(x = karyotype_subtype, y = chr21_index)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, size = 0.7, alpha = 0.5) +
  geom_hline(yintercept = c(1, 1.5), linetype = "dotted") +
  labs(title = "C  chr21 expression index by karyotype subtype", x = NULL, y = "median chr21 ratio to control") + th
panels <- list(pA, pC)
if (!is.null(B)) panels <- c(panels, list(
  ggplot(B, aes(x = reorder(cluster, genes_assoc_t21), y = genes_assoc_t21)) + geom_col(fill = "#2C7BB6") + coord_flip() +
    labs(title = "B  Genes associated with each cell fraction (T21, 5% FDR)", x = NULL, y = "genes") + th))
if (!is.null(D)) panels <- c(panels, list(
  ggplot(D[abs(contribution) > 0.01], aes(x = contribution, y = Gene_name, colour = term)) +
    geom_vline(xintercept = 0, colour = "grey60") + geom_point(size = 1.6) +
    labs(title = "D  Per-covariate attribution of the fold-change shift", x = "log2 contribution", y = NULL) +
    th + theme(legend.position = "bottom", legend.text = element_text(size = 6), legend.title = element_blank()) +
    guides(colour = guide_legend(ncol = 4))))
fig <- wrap_plots(panels, ncol = 2)
h <- 5 * ceiling(length(panels) / 2)
ggsave(run$figure("S1_covariate_evidence.pdf"), fig, width = 12, height = h)
ggsave(run$figure("S1_covariate_evidence.png"), fig, width = 12, height = h, dpi = 200)
writeLines(capture.output(sessionInfo()), run$figure("S1_covariate_evidence_session_info.txt"))
cat(sprintf("  Saved: %s\n=== S1 complete ===\n", run$figure("S1_covariate_evidence.{pdf,png}")))
