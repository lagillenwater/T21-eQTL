# ARCHIVED 2026-09-10: superseded by scripts/10_compare_runs.R and the adjusted run
# (config/runs/adjusted.R). This script approximated composition adjustment on a
# log-CPM scale; the run mechanism does it inside DESeq2. Not expected to run.
# 10_composition_adjusted_neighborhood.R
#
# Purpose: Ask how much of each deviating chr21 gene's T21-vs-Control shift,
#          and of its co-expression neighborhood's shift, remains after the
#          linear association with measured blood cell composition is removed.
#
#          Composition is the HTP CyTOF table (Synapse syn31488783): 20
#          FlowSOM populations as percent of CD45+ CD66low cells, one sample
#          per subject, 372 of the 397 analysis-cohort samples matched by
#          LabID. The CD66low gate excludes granulocytes, so neutrophil
#          fraction is NOT among the covariates; "adjusted" means adjusted for
#          the mononuclear composition that was measured.
#
#          Method. Expression is log2(ploidy-normalized CPM + 1): chr21 counts
#          in T21 samples are divided by 1.5 and the library size excludes
#          chr21, as in script 01, so a chr21 gene's T21-vs-Control difference
#          is its deviation from the trisomy expectation. For every gene the
#          karyotype coefficient is taken from two least-squares fits across
#          the 372 matched samples: expression ~ karyotype (unadjusted) and
#          expression ~ karyotype + the 20 fractions (adjusted). Both are on
#          the same scale, so their ratio is the share of the shift retained
#          when composition is held fixed. This scale compresses shifts for
#          low-expressed genes relative to DESeq2's norm_log2FC (kept beside
#          it for reference); the unadjusted-vs-adjusted comparison is within
#          scale. The partner shift is the median
#          shift of the gene's 20 control-defined partners (the partner sets of
#          script 04, selected in all 95 controls) and is tested against the
#          same correlation-matched null as script 04, rebuilt on each scale.
#          A composition R-squared (within the 282 matched T21 subjects) is
#          drawn beside the neighborhood and cis-variant R-squared of script
#          08, with the background genes of script 09 as reference.
#
#          Everything is descriptive. A shift that shrinks under adjustment is
#          one that co-varies with measured composition; a shift that survives
#          is one that does not. Neither statement identifies a cause.
#
# Inputs:
#   - results/tables/chr21_lane_assignments.csv           (script 04)
#   - results/tables/deviation_variance_explained.csv     (script 08)
#   - results/tables/neighborhood_r2_context.csv          (script 09; background gene list)
#   - data/processed/count_matrix.csv, analysis_cohort.csv
#   - data/synapse/syn31488783/HTP_CyTOF_CD45posCD66low_FlowSOM_cluster_percentage_Synapse.txt
#
# Outputs:
#   - results/tables/composition_adjusted_neighborhood.csv   (per deviating gene)
#   - results/tables/composition_karyotype_difference.csv    (per cluster, T21 vs Control)
#   - results/tables/composition_r2_context.csv              (composition R2, all sets)
#   - results/figures/composition_adjusted_neighborhood.{pdf,png}
#   - results/figures/composition_adjusted_neighborhood_session_info.txt

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
source("scripts/lib/neighborhood.R")
source("scripts/lib/variance_explained.R")
source("scripts/lib/composition_adjust.R")

N_PARTNERS <- 20
N_DRAW     <- 300
SEED       <- 1
PLOIDY     <- 1.5
CYTOF_PATH <- "data/synapse/syn31488783/HTP_CyTOF_CD45posCD66low_FlowSOM_cluster_percentage_Synapse.txt"
OUT_STEM   <- "results/figures/composition_adjusted_neighborhood"

COL_SCALE   <- c(unadjusted = "#C45A00", adjusted = "#762A83")
COL_MEASURE <- c(neighborhood = "#762A83", "cis variant" = "#1B7837", composition = "#2C7BB6")

cat("=== T21-eQTL: Composition-adjusted neighborhood check ===\n\n")

# =============================================================================
# STEP 1: Expression, composition, and the matched sample set
# =============================================================================

cat("Step 1: Loading inputs...\n")
lanes     <- fread("results/tables/chr21_lane_assignments.csv")
s08       <- fread("results/tables/deviation_variance_explained.csv")
s09       <- fread("results/tables/neighborhood_r2_context.csv")
count_mat <- fread("data/processed/count_matrix.csv")
cohort    <- fread("data/processed/analysis_cohort.csv")
if (!file.exists(CYTOF_PATH)) {
  stop("CyTOF table not found; run `bash download_synapse_celltypes.sh` first: ", CYTOF_PATH)
}
cytof <- read_cytof_wide(CYTOF_PATH)

# log2-CPM, partner pool and the 95-control matrix exactly as in script 04, so
# the partner sets and the null are the same ones the lane table carries.
lab_ids <- intersect(setdiff(names(count_mat), c("EnsemblID", "Gene_name", "Chr")),
                     cohort$LabID)
karyotype <- cohort$Karyotype[match(lab_ids, cohort$LabID)]
M <- as.matrix(count_mat[, ..lab_ids])
rownames(M) <- count_mat$Gene_name
is_chr21 <- count_mat$Chr == "chr21"
L <- log2(t(t(M) / colSums(M) * 1e6) + 1)
partner_pool <- rowMeans(M) >= 25 & !is_chr21
L_ctrl <- L[partner_pool, karyotype == "Control"]

matched <- intersect(lab_ids, rownames(cytof))
g_t21   <- karyotype[match(matched, lab_ids)] == "T21"
C       <- cytof[matched, , drop = FALSE]
# Ploidy-normalized scale for the shifts (script 01's normalization, log form).
Y_m     <- ploidy_log2cpm(M[, matched], is_chr21, g_t21, ploidy = PLOIDY)
cat(sprintf("  Cohort %d samples; CyTOF-matched %d (%d T21, %d Control); %d clusters\n",
            length(lab_ids), length(matched), sum(g_t21), sum(!g_t21), ncol(C)))

# =============================================================================
# STEP 2: Composition differs by karyotype (descriptive table)
# =============================================================================

cat("\nStep 2: Cluster fractions by karyotype...\n")
comp_diff <- data.table(
  cluster      = colnames(C),
  mean_t21     = colMeans(C[g_t21, , drop = FALSE]),
  mean_control = colMeans(C[!g_t21, , drop = FALSE]))
comp_diff[, `:=`(diff = mean_t21 - mean_control,
                 log2_ratio = log2(mean_t21 / mean_control))]
comp_diff[, p_wilcox := vapply(colnames(C), function(k)
  wilcox.test(C[g_t21, k], C[!g_t21, k])$p.value, 0)]
comp_diff[, q_bh := p.adjust(p_wilcox, "BH")]
comp_diff <- comp_diff[order(-abs(log2_ratio))]
fwrite(comp_diff, "results/tables/composition_karyotype_difference.csv")
print(comp_diff[, .(cluster, mean_t21 = round(mean_t21, 2), mean_control = round(mean_control, 2),
                    log2_ratio = round(log2_ratio, 2), q_bh = signif(q_bh, 2))])

# =============================================================================
# STEP 3: Gene and partner shifts, unadjusted and composition-adjusted
# =============================================================================

cat("\nStep 3: Shifts on both scales...\n")

kc      <- karyotype_coef(Y_m, g_t21, C)             # both fits contain karyotype
dev_un  <- kc[, "unadjusted"]
dev_adj <- kc[, "adjusted"]

pool_names <- rownames(L_ctrl)
null_un  <- partner_null(L_ctrl, dev_un[pool_names],  n_partners = N_PARTNERS, n_draw = N_DRAW, seed = SEED)
null_adj <- partner_null(L_ctrl, dev_adj[pool_names], n_partners = N_PARTNERS, n_draw = N_DRAW, seed = SEED)
cat(sprintf("  Matched null SD: unadjusted %.3f, adjusted %.3f\n", sd(null_un), sd(null_adj)))

dev_genes <- intersect(lanes[sig_lane %in% c("DE_high", "DE_low"), Gene_name], rownames(L))

one_scale <- function(g, partners, dev, null) {
  gene_lfc    <- unname(dev[g])
  partner_lfc <- median(dev[partners], na.rm = TRUE)
  p <- partner_p(gene_lfc, partner_lfc, null)
  z <- partner_z(partner_lfc, null)
  list(gene_lfc = gene_lfc, partner_lfc = partner_lfc, p_partners = p,
       partner_z = z$z, neighborhood = neighborhood_class(gene_lfc, partner_lfc, p))
}

res <- rbindlist(lapply(dev_genes, function(g) {
  partners <- partner_names(g, L_ctrl, L[g, karyotype == "Control"], N_PARTNERS)
  u <- one_scale(g, partners, dev_un,  null_un)
  a <- one_scale(g, partners, dev_adj, null_adj)
  data.table(
    Gene_name = g,
    gene_dev_unadj = u$gene_lfc, gene_dev_adj = a$gene_lfc,
    gene_dev_retained = a$gene_lfc / u$gene_lfc,
    partner_lfc_unadj = u$partner_lfc, partner_lfc_adj = a$partner_lfc,
    partner_z_unadj = u$partner_z, partner_z_adj = a$partner_z,
    p_partners_unadj = u$p_partners, p_partners_adj = a$p_partners,
    neighborhood_unadj = u$neighborhood, neighborhood_adj = a$neighborhood)
}))

# =============================================================================
# STEP 4: Composition R-squared within T21, with the script 09 reference sets
# =============================================================================

cat("\nStep 4: Composition R-squared within matched T21 subjects...\n")
C_t21 <- C[g_t21, , drop = FALSE]
Y_t21 <- Y_m[, g_t21]
r2_comp_all <- r2_rows_multi(Y_t21, C_t21)
names(r2_comp_all) <- rownames(Y_t21)

ctx_sets <- unique(s09[definition == "observed", .(Gene_name, gene_set)])
ctx_sets <- ctx_sets[Gene_name %in% names(r2_comp_all)]
ctx_sets[, r2_composition := r2_comp_all[Gene_name]]
ctx_sets[, n_t21 := sum(g_t21)]
fwrite(ctx_sets, "results/tables/composition_r2_context.csv")
ctx_summary <- ctx_sets[, .(n = .N, median_r2_composition = round(median(r2_composition, na.rm = TRUE), 3)),
                        by = gene_set]
print(ctx_summary)

res[, r2_composition := r2_comp_all[Gene_name]]
res[, `:=`(n_matched = length(matched), n_matched_t21 = sum(g_t21), n_matched_control = sum(!g_t21))]
res <- merge(res, lanes[, .(Gene_name, sig_lane, tier, eqtl_lane, norm_log2FC,
                            neighborhood_s04 = neighborhood, partner_z_s04 = partner_z)],
             by = "Gene_name")
res <- merge(res, s08[, .(Gene_name, r2_neighborhood, r2_cis_best, r2_cis_gtex_lead)],
             by = "Gene_name", all.x = TRUE)
res <- res[order(sig_lane, -abs(gene_dev_unadj))]
fwrite(res, "results/tables/composition_adjusted_neighborhood.csv")

cat("\nPer-gene shifts (deviation scale, log2) and neighborhood labels:\n")
print(res[, .(Gene_name, sig_lane, dev_unadj = round(gene_dev_unadj, 2), dev_adj = round(gene_dev_adj, 2),
              retained = round(gene_dev_retained, 2), pz_unadj = round(partner_z_unadj, 1),
              pz_adj = round(partner_z_adj, 1), nb_unadj = neighborhood_unadj, nb_adj = neighborhood_adj,
              r2_comp = round(r2_composition, 2), r2_nbhd = round(r2_neighborhood, 2))])
cat(sprintf("\nLabels changed by adjustment: %d of %d\n",
            sum(res$neighborhood_unadj != res$neighborhood_adj), nrow(res)))

# =============================================================================
# STEP 5: Figure
# =============================================================================

cat("\nStep 5: Drawing the figure...\n")

res[, lane := factor(fifelse(sig_lane == "DE_high", "DE high", "DE low"), levels = c("DE high", "DE low"))]
gene_order <- res[order(lane, gene_dev_unadj), Gene_name]
res[, gene := factor(Gene_name, levels = gene_order)]

base_theme <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.5, colour = "grey30"),
        legend.position = "bottom", legend.text = element_text(size = 7),
        legend.key.size = unit(0.75, "lines"))

dumbbell <- function(dt, x_un, x_adj, title, subtitle, xlab, vlines = 0) {
  long <- rbind(dt[, .(gene, lane, scale = "unadjusted", x = get(x_un))],
                dt[, .(gene, lane, scale = "adjusted",   x = get(x_adj))])
  long[, scale := factor(scale, levels = names(COL_SCALE))]
  ggplot(long, aes(x = x, y = gene)) +
    geom_vline(xintercept = vlines, colour = "grey60", linetype = "dotted", linewidth = 0.3) +
    geom_line(aes(group = gene), colour = "grey75", linewidth = 0.5) +
    geom_point(aes(colour = scale), size = 2.2) +
    facet_grid(lane ~ ., scales = "free_y", space = "free_y") +
    scale_colour_manual(values = COL_SCALE, limits = names(COL_SCALE), name = NULL,
                        labels = c(unadjusted = "Unadjusted", adjusted = "Adjusted for 20 CyTOF cluster fractions")) +
    labs(title = title, subtitle = subtitle, x = xlab, y = NULL) +
    base_theme
}

panel_a <- dumbbell(res, "gene_dev_unadj", "gene_dev_adj",
  "A  Gene deviation from the trisomy expectation",
  sprintf("Karyotype coefficient on log2(ploidy-normalized CPM + 1); %d matched samples", length(matched)),
  "Deviation (log2)")

panel_b <- dumbbell(res, "partner_z_unadj", "partner_z_adj",
  "B  Neighborhood shift",
  "Median shift of the 20 control-defined partners, in SDs of the matched null (dotted: 2 SD)",
  "Partner shift (null SD units)", vlines = c(-2, 0, 2))

r2_long <- rbind(
  res[, .(gene, lane, measure = "neighborhood", r2 = r2_neighborhood)],
  res[, .(gene, lane, measure = "cis variant",  r2 = r2_cis_best)],
  res[, .(gene, lane, measure = "composition",  r2 = r2_composition)])
r2_long[, measure := factor(measure, levels = names(COL_MEASURE))]
bg_med <- ctx_summary[gene_set == "background", median_r2_composition]
panel_c <- ggplot(r2_long, aes(x = r2, y = gene, colour = measure)) +
  geom_vline(xintercept = bg_med, colour = COL_MEASURE["composition"], linetype = "dashed", linewidth = 0.35) +
  geom_line(aes(group = gene), colour = "grey85", linewidth = 0.4) +
  geom_point(size = 2.2, na.rm = TRUE) +
  annotate("text", x = bg_med, y = -Inf, vjust = -0.4, hjust = -0.05, size = 2.3,
           colour = COL_MEASURE["composition"],
           label = sprintf("background median\ncomposition R2 = %.2f", bg_med)) +
  facet_grid(lane ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = COL_MEASURE, limits = names(COL_MEASURE), name = NULL,
                      labels = c(neighborhood = "Neighborhood score (script 08)",
                                 "cis variant" = "Best cis variant (script 08)",
                                 composition = "20 CyTOF cluster fractions")) +
  labs(title = "C  Within-T21 R-squared on each measure",
       subtitle = sprintf("Composition: %d matched T21 subjects;\nneighborhood and cis variant: 302 genotyped T21 subjects (script 08)",
                          sum(g_t21)),
       x = expression(R^2 ~ "within T21"), y = NULL) +
  base_theme

fig <- (panel_a | panel_b | panel_c) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(paste0(OUT_STEM, ".pdf"), fig, width = 13, height = 7.5)
ggsave(paste0(OUT_STEM, ".png"), fig, width = 13, height = 7.5, dpi = 200)
cat(sprintf("  Saved: %s.{pdf,png}\n", OUT_STEM))

writeLines(capture.output(sessionInfo()), paste0(OUT_STEM, "_session_info.txt"))
cat("\n=== Composition-adjusted neighborhood check complete ===\n")
