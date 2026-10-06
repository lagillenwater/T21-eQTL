# 19_gtex_tissue_expression.R
#
# Purpose: Are the deviating genes meaningfully expressed in whole blood, and
#          does whole-blood expression track the eQTL results (manuscript
#          issue #30)? A gene that is barely expressed in blood may deviate
#          through noise or a small cell population, and its whole-blood
#          eQTLs may be weak or untested, which bears on the "tested, not
#          detected" calls.
#
#   Every classified chr21 gene (the run's lane table) is placed against the
#   GTEx v10 median TPM of all 68 tissues (scripts/lib/gtex_expression.R):
#   whole-blood TPM, the rank of whole blood among the tissues, the top
#   tissue, tissue specificity (tau), and the percentile of the gene's
#   whole-blood TPM among all classified chr21 genes and among the
#   Expected-dosage genes. GTEx whole blood is not globin-depleted, so it
#   ranks near the bottom of the tissues for most genes; the other chr21
#   genes, measured in the same tissue, are the yardstick, not the rank.
#
#   Against the eQTLs, per gene: whether GTEx tested it in whole blood (its
#   own expression inclusion rule for eQTL mapping, read from the eGenes
#   file), its eGene q and aFC, and the within-T21 gene-level test (script 03
#   for the deviating genes, script 16 for the Expected-dosage genes; each
#   set BH-adjusted on its own). Whole-blood TPM is compared, Wilcoxon
#   rank-sum, between GTEx-tested and untested genes, GTEx eGenes and
#   non-eGenes, and T21-detected and not-detected genes; and a logistic curve
#   of T21 detection on log10 whole-blood TPM, fitted on the tested
#   Expected-dosage genes, gives each tested deviating gene a detection
#   probability at its own expression level. Annotation only: no lane changes.
#
#   The cohort's own expression on the same scale: TPM from the raw counts of
#   the run's T21 participants (exon-union gene lengths from the GENCODE v44
#   GTF), set against GTEx's whole-blood median for every chr21 gene. The two
#   are not on an identical scale: the HTP libraries are globin-depleted and
#   GTEx's are not, so every non-globin gene reads higher in the HTP. The
#   genome-wide offset (non-chr21 genes) is written beside them.
#
# Inputs (per run):
#   - tables/chr21_lane_assignments.csv (script 04), eqtl_gene_level_perm.csv
#     (script 03), base_rate_expected_dosage.csv (script 16; optional)
#   - processed/analysis_cohort.csv (script 01), data/processed/count_matrix.csv
#     (script 00), data/gencode.v44.basic.annotation.gtf.gz
#   - data/gtex/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz
#     (downloaded on first use), data/Whole_Blood.v10.eGenes.txt.gz
# Outputs (per run):
#   - tables/gtex_tissue_expression_genes.csv    one row per classified chr21 gene
#   - tables/gtex_tissue_expression_summary.csv  per lane, and all genes
#   - tables/gtex_tissue_expression_eqtl.csv     whole-blood TPM by eQTL outcome
#   - tables/gtex_tissue_expression_detection.csv  T21 detection against whole-blood TPM
#   - tables/gtex_tissue_expression_scale.csv      HTP T21 TPM against GTEx whole blood, genome-wide
#   - figures/gtex_tissue_expression.{pdf,png}
#
# Usage: Rscript scripts/19_gtex_tissue_expression.R --run adjusted

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})
source("scripts/lib/gtex_expression.R")
source("scripts/lib/run.R"); run <- load_run()

EGENES  <- "data/Whole_Blood.v10.eGenes.txt.gz"
EGENE_Q <- 0.05
# Direction of deviation, as in scripts 07, 11 and 18; every other gene grey.
COL_DIR <- c(`Higher than expected` = "#B2182B", `Lower than expected` = "#2166AC", `Other chr21 gene` = "grey60")
LANE_LABEL <- c(DE_high = "DE high", DE_low = "DE low", Expected_dosage = "Expected dosage",
                Not_DE_outside_noise = "Not DE", Not_assessable = "Not assessable",
                High_repeats = "High repeats", Low_expression = "Low expression")
th <- theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10), panel.grid.minor = element_blank())

cat(sprintf("=== T21-eQTL: GTEx tissue expression of chr21 genes [run: %s] ===\n\n", run$name))

# ---- Tissue expression of every classified chr21 gene -------------------------
genes <- fread(run$table("chr21_lane_assignments.csv"))[
  , .(Gene_name, ensembl_stable, Gene_type, sig_lane, tier, baseMean, eqtl_lane)]
tpm <- read_gtex_median_tpm(gtex_median_tpm_file())
tis <- tissue_expression_summary(tpm[ensembl_stable %in% genes$ensembl_stable])
setnames(tis, sub("^focal_", "wb_", names(tis)))
genes <- merge(genes, tis, by = "ensembl_stable", all.x = TRUE)
genes[, in_gtex := !is.na(wb_tpm)]
pct_in <- function(x, ref) ecdf(ref)(x)
genes[in_gtex == TRUE, `:=`(
  wb_tpm_pct_chr21           = pct_in(wb_tpm, wb_tpm),
  wb_tpm_pct_expected_dosage = pct_in(wb_tpm, wb_tpm[sig_lane == "Expected_dosage"]))]
cat(sprintf("  %d classified chr21 genes, %d in the GTEx v10 median-TPM table (%d tissues)\n",
            nrow(genes), sum(genes$in_gtex), tis$n_tissues[1]))

# ---- The run's T21 participants on the TPM scale ---------------------------------
len <- exon_union_length(read_gtf_exons())
cm  <- fread("data/processed/count_matrix.csv")
cm[, ensembl_stable := sub("[.].*$", "", EnsemblID)]
cm  <- cm[!duplicated(ensembl_stable)]
t21_labs <- intersect(run$cohort()[Karyotype == "T21", LabID], names(cm))
counts <- as.matrix(cm[, ..t21_labs]); rownames(counts) <- cm$ensembl_stable
htp_tpm <- counts_to_tpm(counts, setNames(len$length_bp, len$ensembl_stable))
cat(sprintf("  HTP TPM: %d T21 participants, %d of %d genes with a GENCODE v44 length\n",
            ncol(htp_tpm), nrow(htp_tpm), nrow(counts)))
htp_q <- t(apply(htp_tpm[intersect(genes$ensembl_stable, rownames(htp_tpm)), , drop = FALSE], 1,
                 quantile, probs = c(0.25, 0.5, 0.75), names = FALSE))
genes <- merge(genes, data.table(ensembl_stable = rownames(htp_q), htp_t21_tpm_q25 = htp_q[, 1],
                                 htp_t21_tpm_median = htp_q[, 2], htp_t21_tpm_q75 = htp_q[, 3]),
               by = "ensembl_stable", all.x = TRUE)

# Genome-wide scale check, non-chr21 genes (no trisomy) measured in both:
# HTP T21 median TPM over GTEx whole-blood median TPM, and the share of each
# library taken by the three main hemoglobin genes.
HEMOGLOBIN <- c(HBB = "ENSG00000244734", HBA1 = "ENSG00000206172", HBA2 = "ENSG00000188536")
gw <- data.table(ensembl_stable = rownames(htp_tpm), htp = apply(htp_tpm, 1, median))
gw <- merge(gw, cm[, .(ensembl_stable, Chr)], by = "ensembl_stable")
gw <- merge(gw, tpm[, .(ensembl_stable, gtex = get(GTEX_WHOLE_BLOOD))], by = "ensembl_stable")
gw <- gw[Chr != "chr21" & htp > 0 & gtex > 0]
scale_tbl <- data.table(
  genes_compared = nrow(gw),
  median_htp_over_gtex = round(2^median(log2(gw$htp / gw$gtex)), 2),
  spearman = round(cor(gw$htp, gw$gtex, method = "spearman"), 3),
  hemoglobin_share_htp = signif(median(colSums(htp_tpm[intersect(HEMOGLOBIN, rownames(htp_tpm)), , drop = FALSE])) / 1e6, 3),
  hemoglobin_share_gtex = signif(sum(tpm[ensembl_stable %in% HEMOGLOBIN, get(GTEX_WHOLE_BLOOD)]) / 1e6, 3))
fwrite(scale_tbl, run$table("gtex_tissue_expression_scale.csv"))

# ---- eQTL results: GTEx whole blood and the within-T21 test -------------------
eg <- fread(EGENES)[, .(ensembl_stable = sub("[.].*$", "", gene_id), gtex_qval = qval, gtex_afc = afc)]
genes <- merge(genes, eg, by = "ensembl_stable", all.x = TRUE)
genes[, `:=`(gtex_wb_tested = !is.na(gtex_qval), gtex_egene = !is.na(gtex_qval) & gtex_qval < EGENE_Q)]

t21 <- fread(run$table("eqtl_gene_level_perm.csv"))[
  !is.na(p_gene_perm), .(Gene_name, t21_q_gene_bh = q_gene_bh, t21_detected = cis_eqtl_detected)]
base_rate <- run$table("base_rate_expected_dosage.csv")
if (file.exists(base_rate)) {
  t21 <- rbind(t21, fread(base_rate)[!is.na(p_gene_perm),
                                     .(Gene_name, t21_q_gene_bh = q_gene_bh, t21_detected = detected)])
} else {
  warning("no ", base_rate, " (script 16): Expected-dosage genes carry no T21 test result", call. = FALSE)
}
genes <- merge(genes, t21, by = "Gene_name", all.x = TRUE)
genes[, t21_tested := !is.na(t21_detected)]
genes[, direction := factor(fcase(sig_lane == "DE_high", "Higher than expected",
                                  sig_lane == "DE_low",  "Lower than expected",
                                  default = "Other chr21 gene"), levels = names(COL_DIR))]

# ---- T21 detection against whole-blood TPM ------------------------------------
# Fitted on the tested Expected-dosage genes; read off at each tested
# deviating gene's own whole-blood TPM.
ref <- genes[sig_lane == "Expected_dosage" & t21_tested & in_gtex & wb_tpm > 0]
genes[, p_detect_wb_tpm := NA_real_]
detection_tbl <- NULL
if (nrow(ref) >= 10 && length(unique(ref$t21_detected)) == 2) {
  fit <- glm(t21_detected ~ log10(wb_tpm), family = binomial, data = ref)
  genes[t21_tested & in_gtex & wb_tpm > 0, p_detect_wb_tpm := predict(fit, newdata = .SD, type = "response")]
  d <- genes[sig_lane %in% c("DE_high", "DE_low") & !is.na(p_detect_wb_tpm)]
  detection_tbl <- data.table(
    reference = "Expected dosage, tested in T21", n_reference = nrow(ref),
    detected_reference = sum(ref$t21_detected),
    slope_per_log10_tpm = round(coef(fit)[[2]], 3),
    slope_p = signif(summary(fit)$coefficients[2, 4], 3),
    n_deviating = nrow(d), observed_deviating = sum(d$t21_detected),
    expected_deviating = round(sum(d$p_detect_wb_tpm), 2),
    expected_sd = round(sqrt(sum(d$p_detect_wb_tpm * (1 - d$p_detect_wb_tpm))), 2))
  fwrite(detection_tbl, run$table("gtex_tissue_expression_detection.csv"))
} else {
  warning("too few tested Expected-dosage genes for the detection curve; no detection table", call. = FALSE)
  unlink(run$table("gtex_tissue_expression_detection.csv"))   # never leave a previous run's table behind
}

# ---- Per-gene table -------------------------------------------------------------
lane_order <- names(LANE_LABEL)
genes[, lane_f := factor(sig_lane, levels = lane_order)]
setorder(genes, lane_f, wb_tpm, na.last = TRUE)
fwrite(genes[, .(Gene_name, ensembl_stable, Gene_type, sig_lane, tier, baseMean = round(baseMean, 1), in_gtex,
                 wb_tpm = signif(wb_tpm, 4), wb_rank, n_tissues, top_tissue, top_tpm = signif(top_tpm, 4),
                 wb_frac_of_top = signif(wb_frac_of_top, 3), tau = round(tau, 3),
                 wb_tpm_pct_chr21 = round(wb_tpm_pct_chr21, 3),
                 wb_tpm_pct_expected_dosage = round(wb_tpm_pct_expected_dosage, 3),
                 htp_t21_tpm_median = signif(htp_t21_tpm_median, 4), htp_t21_tpm_q25 = signif(htp_t21_tpm_q25, 4),
                 htp_t21_tpm_q75 = signif(htp_t21_tpm_q75, 4),
                 gtex_wb_tested, gtex_qval = signif(gtex_qval, 3), gtex_egene, gtex_afc = round(gtex_afc, 3),
                 eqtl_lane, t21_tested, t21_q_gene_bh = signif(t21_q_gene_bh, 3), t21_detected,
                 p_detect_wb_tpm = round(p_detect_wb_tpm, 3))],
       run$table("gtex_tissue_expression_genes.csv"))

# ---- Per lane --------------------------------------------------------------------
lane_stats <- function(d) d[, .(
  n = .N, n_in_gtex = sum(in_gtex),
  median_wb_tpm = signif(median(wb_tpm, na.rm = TRUE), 3),
  q25_wb_tpm = signif(quantile(wb_tpm, 0.25, na.rm = TRUE), 3),
  q75_wb_tpm = signif(quantile(wb_tpm, 0.75, na.rm = TRUE), 3),
  median_wb_rank = as.numeric(median(wb_rank, na.rm = TRUE)),
  median_tau = round(median(tau, na.rm = TRUE), 3),
  wb_top_tissue = sum(top_tissue == "Whole_Blood", na.rm = TRUE),
  gtex_wb_tested = sum(gtex_wb_tested), gtex_egenes = sum(gtex_egene),
  t21_tested = sum(t21_tested), t21_detected = sum(t21_detected, na.rm = TRUE))]
lane_tbl <- rbind(genes[, lane_stats(.SD), by = .(sig_lane, lane_f)][order(lane_f)][, lane_f := NULL],
                  genes[, lane_stats(.SD)][, sig_lane := "All classified chr21"], use.names = TRUE)
fwrite(lane_tbl, run$table("gtex_tissue_expression_summary.csv"))

# ---- Whole-blood TPM by eQTL outcome -------------------------------------------
# Each comparison is two groups of genes; Wilcoxon rank-sum on whole-blood TPM.
deviating <- c("DE_high", "DE_low")
groups <- rbind(
  genes[in_gtex == TRUE, .(comparison = "GTEx tested in whole blood (all chr21)",
                           group = ifelse(gtex_wb_tested, "tested", "not tested"), wb_tpm)],
  genes[gtex_wb_tested & in_gtex, .(comparison = "GTEx whole-blood eGene at q < 0.05 (chr21 genes GTEx tested)",
                                    group = ifelse(gtex_egene, "eGene", "not eGene"), wb_tpm)],
  genes[t21_tested & in_gtex, .(comparison = "T21 cis-eQTL detected (Expected dosage + deviating)",
                                group = ifelse(t21_detected, "detected", "not detected"), wb_tpm)],
  genes[t21_tested & in_gtex & sig_lane == "Expected_dosage",
        .(comparison = "T21 cis-eQTL detected (Expected dosage)",
          group = ifelse(t21_detected, "detected", "not detected"), wb_tpm)],
  genes[in_gtex & sig_lane %in% c(deviating, "Expected_dosage"),
        .(comparison = "Deviating vs Expected dosage",
          group = ifelse(sig_lane %in% deviating, "deviating", "Expected dosage"), wb_tpm)])
wilcox_p <- function(x, g) {
  g <- droplevels(as.factor(g))
  if (nlevels(g) == 2 && all(table(g) >= 2)) wilcox.test(x ~ g)$p.value else NA_real_
}
GROUP_ORDER <- c("tested", "not tested", "eGene", "not eGene", "detected", "not detected",
                 "deviating", "Expected dosage")
groups[, `:=`(comparison = factor(comparison, levels = unique(comparison)),
              group = factor(group, levels = GROUP_ORDER))]
eqtl_tbl <- groups[, .(n = .N, median_wb_tpm = signif(median(wb_tpm), 3),
                       q25_wb_tpm = signif(quantile(wb_tpm, 0.25), 3), q75_wb_tpm = signif(quantile(wb_tpm, 0.75), 3)),
                   keyby = .(comparison, group)]
eqtl_tbl <- merge(eqtl_tbl, groups[, .(wilcoxon_p = signif(wilcox_p(wb_tpm, group), 3)), by = comparison],
                  by = "comparison", sort = FALSE)
fwrite(eqtl_tbl, run$table("gtex_tissue_expression_eqtl.csv"))

# ---- Figure ---------------------------------------------------------------------
# Every classified chr21 gene: A whole-blood TPM by lane, B GTEx eQTL strength
# against it (genes tested in T21), C the cohort's TPM against GTEx's. One row
# per deviating gene: D whole blood among the GTEx tissues. Each panel carries
# its own legend.
#
# TPM spans 0 to thousands; a pseudo-log axis keeps the zeros and is
# near-logarithmic from 0.01 TPM up.
tpm_trans  <- scales::pseudo_log_trans(sigma = 0.01, base = 10)
tpm_breaks <- c(0, 0.1, 1, 10, 100, 1000, 10000)
tpm_labels <- c("0", "0.1", "1", "10", "100", "1k", "10k")
tpm_scale   <- scale_x_continuous(trans = tpm_trans, breaks = tpm_breaks, labels = tpm_labels)
tpm_scale_y <- scale_y_continuous(trans = tpm_trans, breaks = tpm_breaks, labels = tpm_labels)
legend_bottom <- theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0))
pd <- genes[in_gtex == TRUE]
n_lane <- pd[, .N, by = sig_lane]
pd[, lane_lab := factor(sprintf("%s (n = %d)", LANE_LABEL[sig_lane], n_lane[match(pd$sig_lane, sig_lane), N]),
                        levels = rev(sprintf("%s (n = %d)", LANE_LABEL[lane_order],
                                             n_lane[match(lane_order, sig_lane), N])))]
pd[, lane_lab := droplevels(lane_lab)]
pd[, label := ifelse(sig_lane %in% deviating, Gene_name, "")]
jit <- position_jitter(width = 0, height = 0.18, seed = 19)
pA <- ggplot(pd, aes(wb_tpm, lane_lab)) +
  geom_boxplot(outlier.shape = NA, fill = "grey95", colour = "grey55", linewidth = 0.35, width = 0.6) +
  geom_point(data = pd[!sig_lane %in% deviating], position = jit, colour = COL_DIR[["Other chr21 gene"]],
             size = 1, alpha = 0.7) +
  geom_point(data = pd[sig_lane %in% deviating], aes(colour = direction), position = jit, size = 2.4) +
  geom_text_repel(data = pd[sig_lane %in% deviating], aes(label = label), position = jit, size = 2.4,
                  colour = "grey15", min.segment.length = 0, segment.size = 0.2, max.overlaps = 40) +
  tpm_scale +
  scale_colour_manual(values = COL_DIR[1:2], name = "Deviating gene") +
  labs(x = "GTEx v10 whole-blood median TPM", y = NULL,
       title = "A  Whole-blood expression of chr21 genes, by lane",
       subtitle = "Classified chr21 genes in GTEx v10; box: median and IQR") +
  th + legend_bottom

# Genes tested in T21 (scripts 03, 16); all of them are GTEx-tested, since the
# T21 test runs on GTEx cis variants.
pb <- pd[t21_tested == TRUE & gtex_wb_tested == TRUE]
pb[, t21 := factor(ifelse(t21_detected, "detected", "not detected"), levels = c("detected", "not detected"))]
pB <- ggplot(pb, aes(wb_tpm, -log10(gtex_qval))) +
  geom_hline(yintercept = -log10(EGENE_Q), linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_point(data = pb[!sig_lane %in% deviating], aes(shape = t21), colour = COL_DIR[["Other chr21 gene"]],
             size = 1.4, stroke = 0.6) +
  geom_point(data = pb[sig_lane %in% deviating], aes(shape = t21, colour = direction), size = 2.4, stroke = 0.9) +
  geom_text_repel(data = pb[sig_lane %in% deviating], aes(label = label), size = 2.4, colour = "grey15",
                  min.segment.length = 0, segment.size = 0.2, max.overlaps = 40) +
  tpm_scale +
  scale_colour_manual(values = COL_DIR[1:2], name = "Deviating gene") +
  scale_shape_manual(values = c(detected = 16, `not detected` = 1), name = "cis-eQTL in T21") +
  labs(x = "GTEx v10 whole-blood median TPM", y = expression(-log[10] ~ "GTEx whole-blood eGene q"),
       title = "B  GTEx whole-blood eQTL against expression",
       subtitle = sprintf("chr21 genes tested in T21 (n = %d); grey: Expected dosage; dashed: q = %.2f",
                          nrow(pb), EGENE_Q)) +
  th + legend_bottom

# C: the cohort against GTEx for every classified chr21 gene. The HTP libraries
# are globin-depleted, GTEx's are not, so the genome-wide offset (non-chr21
# genes, scale_tbl) is drawn beside the line of equality.
sc <- pd[!is.na(htp_t21_tpm_median)]
off <- scale_tbl$median_htp_over_gtex
ref_x <- c(0, 10^seq(-2, 4, length.out = 200))
ref_lines <- rbind(data.table(x = ref_x, y = ref_x, line = "equal"),
                   data.table(x = ref_x, y = off * ref_x, line = sprintf("genome-wide HTP/GTEx, %.1fx", off)))
ref_lines[, line := factor(line, levels = unique(line))]
sc_max <- max(sc$wb_tpm, sc$htp_t21_tpm_median)
pC <- ggplot(sc, aes(wb_tpm, htp_t21_tpm_median)) +
  geom_line(data = ref_lines[x <= sc_max & y <= sc_max], aes(x, y, linetype = line), colour = "grey40",
            linewidth = 0.4) +
  geom_point(data = sc[!sig_lane %in% deviating], colour = COL_DIR[["Other chr21 gene"]], size = 1, alpha = 0.7) +
  geom_point(data = sc[sig_lane %in% deviating], aes(colour = direction), size = 2.4) +
  geom_text_repel(data = sc[sig_lane %in% deviating], aes(label = label), size = 2.4, colour = "grey15",
                  min.segment.length = 0, segment.size = 0.2, max.overlaps = 40) +
  tpm_scale + tpm_scale_y +
  scale_colour_manual(values = COL_DIR[1:2], name = "Deviating gene") +
  scale_linetype_manual(values = c("dotted", "dashed"), name = NULL) +
  labs(x = "GTEx v10 whole-blood median TPM", y = sprintf("T21 participants' median TPM (n = %d)", ncol(htp_tpm)),
       title = "C  The cohort against GTEx whole blood",
       subtitle = sprintf("Classified chr21 genes (n = %d), Spearman %.2f; grey: other chr21 genes",
                          nrow(sc), cor(sc$wb_tpm, sc$htp_t21_tpm_median, method = "spearman"))) +
  th + legend_bottom

# D: one row per deviating gene, split by direction.
dev <- genes[sig_lane %in% deviating & in_gtex == TRUE]
setorder(dev, sig_lane, wb_tpm)
dev[, `:=`(gene_f = factor(Gene_name, levels = Gene_name), gene_i = seq_len(.N),
           dir_f  = factor(LANE_LABEL[sig_lane], levels = LANE_LABEL[deviating]))]
tissues <- setdiff(names(tpm), c("ensembl_stable", "gene_name"))
pt <- melt(tpm[ensembl_stable %in% dev$ensembl_stable], id.vars = "ensembl_stable", measure.vars = tissues,
           variable.name = "tissue", value.name = "tpm")
pt <- merge(pt, dev[, .(ensembl_stable, gene_i, dir_f)], by = "ensembl_stable")
pt[, group := gtex_tissue_group(as.character(tissue))]
# Each tissue group gets its own sub-lane within the gene's row, whole blood
# in the middle (on the gene's label), three groups above and three below, so
# position separates the groups and colour only has to tell neighbouring
# sub-lanes apart (adjacent-pair validated palette).
COL_TISSUE <- c(`Whole blood` = "grey10", `Immune and lymphoid` = "#e87ba4", `Brain and nerve` = "#4a3aa7",
                `Heart, vessels and muscle` = "#eb6834", Digestive = "#1baf7a", Reproductive = "#008300",
                Other = "#8c8c8c")
sublane_order <- c("Immune and lymphoid", "Brain and nerve", "Heart, vessels and muscle", "Whole blood",
                   "Digestive", "Reproductive", "Other")
lane_offset <- setNames(seq(0.33, -0.33, length.out = length(sublane_order)), sublane_order)
pt[, y := gene_i + lane_offset[as.character(group)]]
rank_lab <- sprintf("%s (%s/%d)", dev$Gene_name, format(dev$wb_rank, trim = TRUE), dev$n_tissues)
n_in_group <- table(gtex_tissue_group(tissues))
group_lab <- setNames(sprintf("%s (%d)", names(n_in_group), as.integer(n_in_group)), names(n_in_group))
x_max <- max(pt$tpm)
gene_rows <- scale_y_continuous(breaks = dev$gene_i, labels = rank_lab, expand = expansion(add = 0.17))
pD <- ggplot(pt, aes(tpm, y, colour = group)) +
  geom_hline(data = dev[gene_i > min(gene_i), .(dir_f, yint = gene_i - 0.5), by = dir_f][, -1],
             aes(yintercept = yint), colour = "grey85", linewidth = 0.3) +
  geom_point(data = pt[group != "Whole blood"], size = 1.1, alpha = 0.85) +
  geom_point(data = pt[group == "Whole blood"], size = 2.8) +
  facet_grid(dir_f ~ ., scales = "free_y", space = "free_y") +
  tpm_scale + gene_rows + coord_cartesian(xlim = c(0, x_max)) +
  scale_colour_manual(values = COL_TISSUE, labels = group_lab, name = "GTEx tissue (n)", drop = FALSE) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = c(2.8, rep(2, 6)), alpha = 1))) +
  labs(x = "GTEx v10 median TPM", y = NULL,
       title = "D  Whole blood among the GTEx tissues, deviating genes",
       subtitle = sprintf(paste0("One point per tissue; whole blood (large) on the gene's line, other tissues in a\n",
                                 "sub-lane per group; label: rank of whole blood among the %d tissues"),
                          length(tissues))) +
  th + legend_bottom + theme(strip.text = element_text(size = 8), panel.grid.major.y = element_blank())

fig <- (pA + pB) / (pC + pD) + plot_layout(heights = c(1, 1.25))
stem <- run$figure("gtex_tissue_expression")
ggsave(paste0(stem, ".pdf"), fig, width = 14, height = 13)
ggsave(paste0(stem, ".png"), fig, width = 14, height = 13, dpi = 300)

# ---- Console -------------------------------------------------------------------
cat("\nPer lane:\n"); print(lane_tbl)
cat("\nWhole-blood TPM by eQTL outcome:\n"); print(eqtl_tbl)
if (!is.null(detection_tbl)) { cat("\nT21 detection against whole-blood TPM:\n"); print(detection_tbl) }
cat("\nHTP T21 TPM against GTEx whole blood, genome-wide:\n"); print(scale_tbl)
cat("\nDeviating genes:\n")
print(genes[sig_lane %in% deviating, .(Gene_name, sig_lane, baseMean = round(baseMean), wb_tpm = signif(wb_tpm, 3),
                                       wb_rank, top_tissue, tau = round(tau, 2),
                                       pct_chr21 = round(wb_tpm_pct_chr21, 2),
                                       pct_expected = round(wb_tpm_pct_expected_dosage, 2),
                                       htp_t21 = signif(htp_t21_tpm_median, 3),
                                       gtex_wb_tested, gtex_egene, eqtl_lane, p_detect_wb_tpm = round(p_detect_wb_tpm, 2))])
writeLines(capture.output(sessionInfo()), run$table("gtex_tissue_expression_session_info.txt"))
cat("\nDone.\n")
