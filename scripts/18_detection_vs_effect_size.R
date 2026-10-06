# 18_detection_vs_effect_size.R
#
# Purpose: Is a detected cis-eQTL more, or less, common among the deviating
#          genes than among chr21 genes at the expected dosage, once the size
#          of the eQTL is accounted for? Detection by the gene-level test
#          (script 03, and script 16 for the Expected-dosage genes: the same
#          test, variant universe, subjects, artifact and thresholds) depends
#          strongly on effect size, and the deviating genes carry larger GTEx
#          effects than the typical Expected-dosage gene, so their raw
#          detection rates are not comparable.
#
#   The positive set is every Expected-dosage gene that is a GTEx
#   whole-blood eGene (q < 0.05) and was tested. A logistic curve of
#   detection on log2 GTEx |aFC| (lead variant) is fitted to that set; each
#   deviating gene's detection probability is read off it at the gene's own
#   GTEx |aFC|, and the observed number of detected deviating genes is
#   compared with the Poisson-binomial distribution of those probabilities.
#   Sensitivity: the same with log10 baseMean in the curve. A joint model
#   with a deviating indicator gives the odds ratio at matched effect size.
#
#   Effect size, not only detection: each dataset's own top variant. The T21
#   allelic fold change at the gene's T21 best variant (the variant the test
#   selected; ploidy-3 curve, scripts/lib/allelic_fc.R, run covariates fitted
#   jointly) is divided by GTEx's |aFC| at its lead variant (the variant GTEx
#   selected, and the only one where GTEx reports an aFC), so both effects
#   carry their own dataset's selection. The two are often different variants
#   tagging one signal, so the T21 aFC is signed by GTEx's own direction at
#   the T21 best variant: 1 = the GTEx effect, 0 = none, negative = against
#   GTEx. Deviating genes are compared with the Expected-dosage eGenes whose
#   GTEx |aFC| falls in the deviating genes' range (matched effect size),
#   Wilcoxon rank-sum. The ratio with T21 also at the GTEx lead variant is
#   kept as a column (ratio_lead).
#
#   The figure: A, that ratio for the two sets; B, per deviating gene, the
#   spike-in power of the test at its own GTEx aFC (script 15).
#
# Inputs (per run):
#   - tables/base_rate_expected_dosage.csv (script 16), eqtl_gene_level_perm.csv
#     (script 03), chr21_lane_assignments.csv (script 04), spike_in_power.csv
#     and spike_in_power_genes.csv (script 15)
#   - processed/base_rate_genotypes.csv (script 16), genotypes_filtered.csv
#     (script 02), expression_adjusted.csv, analysis_cohort.csv (script 01)
#   - data/Whole_Blood.v10.eGenes.txt.gz, the GTEx allpairs chr21 parquet
# Outputs (per run):
#   - tables/detection_vs_effect_size.csv          one row per tested gene
#   - tables/detection_vs_effect_size_summary.csv  detection: observed vs expected, odds ratio
#   - tables/top_variant_effect_ratio.csv          T21 top-variant / GTEx lead aFC, per gene
#   - tables/top_variant_effect_ratio_summary.csv  that ratio by set, matched effect size
#   - figures/eqtl_detection_overview.{pdf,png}
#
# Usage: Rscript scripts/18_detection_vs_effect_size.R --run adjusted

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(arrow)
})
for (f in c("cohort", "eqtl_fit", "allelic_fc", "covariates", "alleles")) source(sprintf("scripts/lib/%s.R", f))
source("scripts/lib/run.R"); run <- load_run()
source("scripts/lib/figure_style.R")   # strip_titles: panel letters only
TITLES <- file.path(run$dirs$figures, "removed_titles.txt")

EGENES  <- "data/Whole_Blood.v10.eGenes.txt.gz"
ALLPAIRS <- paste0("data/GTEx_Analysis_v10_QTLs_GTEx_Analysis_v10_eQTL",
                   "_all_associations_Whole_Blood.v10.allpairs.chr21.parquet")
EGENE_Q <- 0.05
# Direction of deviation, as in scripts 07 and 11; Expected-dosage genes grey.
COL_DIR <- c(`Higher than expected` = "#B2182B", `Lower than expected` = "#2166AC", `Expected dosage` = "grey55")
th <- theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10), panel.grid.minor = element_blank())

cat(sprintf("=== T21-eQTL: detection against eQTL effect size [run: %s] ===\n\n", run$name))

# ---- Tested genes and their GTEx effect sizes --------------------------------
lanes <- fread(run$table("chr21_lane_assignments.csv"))[, .(Gene_name, ensembl_stable, baseMean, sig_lane)]
eg <- fread(EGENES)[, .(ensembl_stable = sub("\\..*$", "", gene_id), gtex_qval = qval, gtex_afc = afc,
                         lead_variant = variant_id)]
expd <- fread(run$table("base_rate_expected_dosage.csv"))[!is.na(p_gene_perm),
                                                          .(Gene_name, n_variants, q_gene_bh, detected, best_variant)]
dev  <- fread(run$table("eqtl_gene_level_perm.csv"))[!is.na(p_gene_perm),
                                                     .(Gene_name, n_variants, q_gene_bh, detected = cis_eqtl_detected,
                                                       best_variant)]
genes <- rbind(expd[, set := "Expected dosage"], dev[, set := "Deviating"])
genes <- merge(merge(genes, lanes, by = "Gene_name"), eg, by = "ensembl_stable", all.x = TRUE)
genes[, `:=`(egene = !is.na(gtex_qval) & gtex_qval < EGENE_Q, abs_afc = abs(gtex_afc))]
cat(sprintf("  Tested: %d Expected-dosage genes (%d GTEx eGenes), %d deviating genes (%d GTEx eGenes)\n",
            genes[set == "Expected dosage", .N], genes[set == "Expected dosage", sum(egene)],
            genes[set == "Deviating", .N], genes[set == "Deviating", sum(egene)]))
genes <- genes[egene == TRUE]

# ---- Detection curve on the Expected-dosage eGenes ----------------------------
ref <- genes[set == "Expected dosage"]
logit_afc  <- glm(detected ~ log2(abs_afc), family = binomial, data = ref)
logit_both <- glm(detected ~ log2(abs_afc) + log10(baseMean), family = binomial, data = ref)
genes[, `:=`(p_afc  = predict(logit_afc,  newdata = genes, type = "response"),
             p_both = predict(logit_both, newdata = genes, type = "response"))]

#' Distribution of the number of successes among independent Bernoulli(p_i).
poisson_binomial <- function(p) Reduce(function(d, pi) c(d * (1 - pi), 0) + c(0, d * pi), p, 1)
compare <- function(label, p, observed) {
  d <- poisson_binomial(p)
  k <- seq_along(d) - 1L
  data.table(model = label, n_deviating = length(p), observed = observed,
             expected = round(sum(p), 2), expected_sd = round(sqrt(sum(p * (1 - p))), 2),
             p_at_least = signif(sum(d[k >= observed]), 3), p_at_most = signif(sum(d[k <= observed]), 3))
}
d <- genes[set == "Deviating"]
summary_tbl <- rbind(compare("GTEx |aFC|", d$p_afc, sum(d$detected)),
                     compare("GTEx |aFC| + baseMean", d$p_both, sum(d$detected)))

# Odds ratio of detection for deviating vs Expected-dosage at matched effect size
joint <- glm(detected ~ log2(abs_afc) + deviating, family = binomial,
             data = genes[, .(detected, abs_afc, deviating = set == "Deviating")])
ci <- suppressMessages(confint(joint))["deviatingTRUE", ]
summary_tbl[, `:=`(odds_ratio_deviating = c(round(exp(coef(joint)[["deviatingTRUE"]]), 2), NA),
                   or_ci_lo = c(round(exp(ci[[1]]), 2), NA), or_ci_hi = c(round(exp(ci[[2]]), 2), NA))]
summary_tbl[, `:=`(n_expected_dosage_egenes = nrow(ref), detected_expected_dosage = sum(ref$detected),
                   afc_slope_per_doubling = c(round(coef(logit_afc)[[2]], 2), round(coef(logit_both)[[2]], 2)))]
fwrite(summary_tbl, run$table("detection_vs_effect_size_summary.csv"))

bins <- quantile(ref$abs_afc, 0:4 / 4)
ref[, afc_bin := cut(abs_afc, bins, include.lowest = TRUE)]
binned <- ref[, .(n = .N, detected = sum(detected), mid = exp(mean(log(abs_afc)))), by = afc_bin][order(afc_bin)]
binned[, rate := detected / n]
wilson <- function(x, n, z = 1.96) {
  p <- x / n; den <- 1 + z^2 / n; c0 <- (p + z^2 / (2 * n)) / den
  hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  list(c0 - hw, c0 + hw)
}
binned[, c("lo", "hi") := wilson(detected, n)]
out <- genes[, .(Gene_name, set, gtex_afc = round(gtex_afc, 3), gtex_qval = signif(gtex_qval, 3),
                 baseMean = round(baseMean, 1), n_variants, q_gene_bh = signif(q_gene_bh, 3), detected,
                 p_detect_afc = round(p_afc, 3), p_detect_afc_basemean = round(p_both, 3))]
fwrite(out[order(set, -abs(gtex_afc))], run$table("detection_vs_effect_size.csv"))

# ---- T21 allelic fold change ---------------------------------------------------
t21 <- run$cohort()[Karyotype == "T21"]
t21[, subject_id := subject_id_from_labid(LabID)]
COV <- run_covariate_matrix(run, t21)
if (!is.null(COV)) rownames(COV) <- t21$subject_id
E <- run$expression()
geno <- rbind(fread(run$processed("base_rate_genotypes.csv")),
              fread(run$processed("genotypes_filtered.csv"))[karyotype == "T21", .(variant_id, subject_id, alt_dosage)])
geno <- geno[!is.na(alt_dosage) & subject_id %in% t21$subject_id,
             .(alt_dosage = mean(alt_dosage)), by = .(variant_id, subject_id)]
lab_of <- setNames(t21$LabID, t21$subject_id)
t21_afc_at <- function(gene, variant) {
  d <- geno[variant_id == variant]
  if (!nrow(d) || !gene %in% rownames(E)) return(list(t21_afc = NA_real_, t21_afc_se = NA_real_))
  f <- fit_afc(d$alt_dosage, E[gene, lab_of[d$subject_id]], T21_CHR21_PLOIDY,
               covariates = if (is.null(COV)) NULL else COV[d$subject_id, , drop = FALSE])
  list(t21_afc = f$afc, t21_afc_se = f$se)
}
genes[, c("t21_afc", "t21_afc_se") := rbindlist(Map(t21_afc_at, Gene_name, lead_variant))]
genes[, ratio_lead := t21_afc / gtex_afc]

# Each dataset's own top variant: the T21 best variant (the one the test
# selected) against GTEx's lead (the one GTEx selected), so both effects carry
# their own dataset's selection. They are often different variants tagging one
# signal, possibly in opposite phase, so the T21 aFC is signed by GTEx's own
# direction at the T21 best variant (allpairs slope): positive = T21 agrees
# with GTEx. The ratio is then over GTEx's |aFC| at its lead.
genes[, c("t21_afc_best", "t21_afc_best_se") := rbindlist(Map(t21_afc_at, Gene_name, best_variant))]
ap <- as.data.table(read_parquet(ALLPAIRS, col_select = c("gene_id", "variant_id", "slope")))
ap[, ensembl_stable := sub("\\..*$", "", gene_id)]
genes <- merge(genes, ap[, .(ensembl_stable, best_variant = variant_id, gtex_slope_at_best = slope)],
               by = c("ensembl_stable", "best_variant"), all.x = TRUE)
genes[, ratio := sign(gtex_slope_at_best) * t21_afc_best / abs(gtex_afc)]
afc_range <- genes[set == "Deviating", range(abs_afc)]
genes[, matched := set == "Deviating" | (abs_afc >= afc_range[1] & abs_afc <= afc_range[2])]
cmp <- genes[matched == TRUE & !is.na(ratio)]
ratio_tbl <- cmp[, .(n = .N, median_ratio = round(median(ratio), 2), q25 = round(quantile(ratio, 0.25), 2),
                     q75 = round(quantile(ratio, 0.75), 2), detected = sum(detected)), by = set]
ratio_tbl[, `:=`(gtex_abs_afc_range = sprintf("%.2f-%.2f", afc_range[1], afc_range[2]),
                 wilcoxon_p = signif(wilcox.test(ratio ~ set, data = cmp)$p.value, 2))]
ratio_tbl[, median_ratio_lead_only := cmp[!is.na(ratio_lead), round(median(ratio_lead), 2), by = set][match(ratio_tbl$set, set), V1]]
fwrite(ratio_tbl, run$table("top_variant_effect_ratio_summary.csv"))
fwrite(genes[, .(Gene_name, set, lead_variant, best_variant, best_is_lead = best_variant == lead_variant,
                 gtex_afc_lead = round(gtex_afc, 3), gtex_slope_at_best = signif(gtex_slope_at_best, 3),
                 t21_afc_best = round(t21_afc_best, 3), t21_afc_best_se = round(t21_afc_best_se, 3),
                 ratio_top = round(ratio, 3), t21_afc_lead = round(t21_afc, 3), ratio_lead = round(ratio_lead, 3),
                 matched, detected)][order(set, -ratio_top)],
       run$table("top_variant_effect_ratio.csv"))

# ---- Figure ------------------------------------------------------------------
cmp[, set_f := factor(set, levels = c("Deviating", "Expected dosage"),
                      labels = c(sprintf("Deviating\n(n = %d)", sum(cmp$set == "Deviating")),
                                 sprintf("Expected dosage, matched\n(n = %d)", sum(cmp$set == "Expected dosage"))))]
cmp[, outcome := factor(ifelse(detected, "cis-eQTL detected", "tested, not detected"),
                         levels = c("cis-eQTL detected", "tested, not detected"))]
cmp[, direction := factor(fcase(sig_lane == "DE_high", "Higher than expected", sig_lane == "DE_low", "Lower than expected",
                                default = "Expected dosage"), levels = names(COL_DIR))]
cmp[, label := ifelse(set == "Deviating", Gene_name, "")]
jit <- position_jitter(width = 0.12, height = 0, seed = 18)
pA <- ggplot(cmp, aes(set_f, ratio)) +
  geom_hline(yintercept = c(0, 1), linetype = c("dotted", "dashed"), colour = "grey50", linewidth = 0.4) +
  geom_violin(fill = "grey93", colour = "grey60", linewidth = 0.35, scale = "width", width = 0.75, trim = TRUE) +
  stat_summary(fun = median, geom = "crossbar", width = 0.35, linewidth = 0.3, colour = "grey30") +
  geom_point(aes(colour = direction, shape = outcome), position = jit, size = 2.2, stroke = 0.8, fill = "white") +
  geom_text_repel(aes(label = label), position = jit, size = 2.4, colour = "grey15", fontface = "italic",
                  min.segment.length = 0,
                  segment.size = 0.2, max.overlaps = 30) +
  scale_colour_manual(values = COL_DIR, name = NULL, drop = FALSE) +
  scale_shape_manual(values = c(`cis-eQTL detected` = 16, `tested, not detected` = 1), name = "Within T21") +
  labs(x = NULL, y = "T21 / GTEx log2 aFC, each at its own top variant",
       title = "A  Top cis-eQTL effect, T21 against GTEx",
       subtitle = sprintf(paste0("Expected-dosage eGenes matched on GTEx |log2 aFC| %.2f-%.2f\n",
                                 "Median %.2f vs %.2f, Wilcoxon p = %.2g; dashed: equal; dotted: none"),
                          afc_range[1], afc_range[2], ratio_tbl[set == "Deviating", median_ratio],
                          ratio_tbl[set == "Expected dosage", median_ratio], ratio_tbl$wilcoxon_p[1])) +
  th + theme(legend.position = "right")

dir_of <- setNames(as.character(cmp[set == "Deviating", direction]), cmp[set == "Deviating", Gene_name])
sp  <- fread(run$table("spike_in_power.csv"))[, direction := factor(dir_of[Gene_name], levels = names(COL_DIR))]
spg <- fread(run$table("spike_in_power_genes.csv"))[, .(Gene_name, gtex_afc = abs(gtex_afc), power_at_gtex_afc_p05,
                                                        detected_observed)]
spg[, `:=`(direction = factor(dir_of[Gene_name], levels = names(COL_DIR)),
           outcome = factor(ifelse(detected_observed, "cis-eQTL detected", "tested, not detected"),
                            levels = c("cis-eQTL detected", "tested, not detected")))]
pB <- ggplot(sp, aes(afc, power_p05, group = Gene_name, colour = direction)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_line(alpha = 0.6, linewidth = 0.5) +
  geom_point(data = spg, mapping = aes(gtex_afc, power_at_gtex_afc_p05, colour = direction, shape = outcome),
             inherit.aes = FALSE, size = 2.3, stroke = 0.8, fill = "white") +
  geom_text_repel(data = spg, aes(gtex_afc, power_at_gtex_afc_p05, label = Gene_name), inherit.aes = FALSE,
                  fontface = "italic",
                  size = 2.4, colour = "grey15", min.segment.length = 0, segment.size = 0.2, max.overlaps = 30) +
  scale_colour_manual(values = COL_DIR[1:2], name = NULL) +
  scale_shape_manual(values = c(`cis-eQTL detected` = 16, `tested, not detected` = 1), name = "Within T21") +
  scale_x_continuous(limits = c(0, max(sp$afc))) +
  labs(x = "Spiked allelic effect, |log2 aFC|", y = "Power (gene-level p < 0.05)",
       title = "B  Power of the test for each deviating gene",
       subtitle = "Spike-in on the gene's own genotypes and noise; points: at its GTEx |log2 aFC|") +
  th + theme(legend.position = "right")

pA <- strip_titles(pA, "eqtl_detection_overview", "A", TITLES, tag = TRUE)
pB <- strip_titles(pB, "eqtl_detection_overview", "B", TITLES, tag = TRUE)
fig <- pA + pB + plot_layout(widths = c(1, 1.15))
stem <- run$figure("eqtl_detection_overview")
ggsave(paste0(stem, ".pdf"), fig, width = 11, height = 5)
ggsave(paste0(stem, ".png"), fig, width = 11, height = 5, dpi = 300)

cat("\nT21 top-variant aFC / GTEx lead-variant aFC, matched effect size:\n"); print(ratio_tbl)
cat("\nDetection rate of the Expected-dosage GTEx eGenes by GTEx |aFC| quartile:\n"); print(binned[, .(afc_bin, n, detected, rate = round(rate, 2), lo = round(lo, 2), hi = round(hi, 2))])
cat("\nDeviating genes: observed against expected at matched effect size:\n"); print(summary_tbl)
cat("\nPer deviating gene:\n"); print(out[set == "Deviating", .(Gene_name, gtex_afc, detected, p_detect_afc, p_detect_afc_basemean)])
writeLines(capture.output(sessionInfo()), run$table("detection_vs_effect_size_session_info.txt"))
cat("\nDone.\n")
