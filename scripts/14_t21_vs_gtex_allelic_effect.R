# 14_t21_vs_gtex_allelic_effect.R
#
# Purpose: Does an allele have the same per-chromosome effect on a trisomic
#          background as in euploid blood? (manuscript review, concern 1
#          follow-up.) The allele-frequency argument that a cis-eQTL cannot
#          move the T21/D21 ratio assumes it does; if the effect were
#          stronger or weaker in T21, the ratio could move even at equal
#          allele frequencies.
#
#   For every dosage panel of script 11 (deviating genes and positive controls
#   on their own best variant of the within-T21 test; negative controls on the
#   best variant of their decoy set) this script estimates the allelic fold
#   change (aFC, log2 expression of a chromosome carrying the plotted allele
#   over one carrying the other) in T21 and sets it against GTEx whole blood.
#
#   Direction: every effect is per copy of the allele on that panel's x-axis
#   in the dosage boxplots (the deviation-matching allele for deviating genes,
#   the minor allele for controls; scripts/lib/eqtl_figures.R::orient_panels),
#   so a positive aFC is a panel that trends up.
#
#   T21: the aFC curve fitted with ploidy 3 (scripts/lib/allelic_fc.R), not a
#        straight line on dosage, so it is the same quantity GTEx reports at
#        ploidy 2.
#   GTEx: aFC is published only at each gene's lead variant. At the best
#        variant it is the GTEx normalised slope converted with the lead
#        variant's own aFC-to-slope ratio (exact when best = lead). The lead
#        variant itself, where no conversion is needed, is fitted in T21 too
#        and written as a sensitivity table.
#   Negative controls: a decoy variant is >= 5 Mb from the gene, outside the
#        GTEx cis window, so GTEx has no estimate; the expectation is 0.
#
# Inputs (per run):
#   - tables/chr21_lane_assignments.csv, eqtl_gene_level_perm.csv,
#     eqtl_control_negative.csv, eqtl_control_positive.csv
#   - processed/eqtl_target_variants.csv (alleles, GTEx slopes),
#     genotypes_filtered.csv, expression_adjusted.csv, analysis_cohort.csv
#   - data/Whole_Blood.v10.eGenes.txt.gz (GTEx lead variant, afc, afc_se)
# Outputs (per run):
#   - tables/t21_vs_gtex_afc.csv               one row per panel (+ untestable genes)
#   - tables/t21_vs_gtex_afc_lead_variant.csv  same comparison at the GTEx lead variant
#   - tables/t21_vs_gtex_afc_summary.csv       agreement per panel group
#   - figures/t21_vs_gtex_allelic_effect.{pdf,png}  A/B, plus C (spike-in power,
#     script 15) and D (base rate, script 16) when their tables exist
#
# Usage: Rscript scripts/14_t21_vs_gtex_allelic_effect.R --run adjusted

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})
source("scripts/lib/cohort.R")        # subject_id_from_labid
source("scripts/lib/eqtl_fit.R")      # expr_from_matrix
source("scripts/lib/alleles.R")       # minor_dosage, align_slope_to_minor, T21_CHR21_PLOIDY
source("scripts/lib/eqtl_figures.R")  # panel_variants, orient_panels, panel_dosage
source("scripts/lib/allelic_fc.R")
source("scripts/lib/run.R"); run <- load_run()

EGENES  <- "data/Whole_Blood.v10.eGenes.txt.gz"
Z95     <- qnorm(0.975)
ALPHA   <- 0.05
SET_LAB <- c("DE high", "DE low", "Positive control", "Negative control")
COL_SET <- setNames(c("#B2182B", "#2166AC", "#1B7837", "grey50"), SET_LAB)
SHP_SET <- setNames(c(16, 16, 17, 18), SET_LAB)
COL_DIR <- c(up = "#B2182B", down = "#2166AC")

cat(sprintf("=== T21-eQTL: allelic effect in T21 vs GTEx [run: %s] ===\n\n", run$name))

# =============================================================================
# STEP 1: Panels and the allele each is drawn on (as in script 11)
# =============================================================================

lanes <- fread(run$table("chr21_lane_assignments.csv"))
perm  <- fread(run$table("eqtl_gene_level_perm.csv"))
neg   <- fread(run$table("eqtl_control_negative.csv"))
pos   <- fread(run$table("eqtl_control_positive.csv"))
tv    <- fread(run$processed("eqtl_target_variants.csv"))

panels <- panel_variants(perm, lanes, neg, pos)
tv[, gtex_slope_minor := align_slope_to_minor(slope, alt_is_minor)]
variant_alleles <- unique(tv[, .(variant_id, REF, ALT, minor_allele, major_allele,
                                 alt_is_minor, gtex_maf)], by = "variant_id")
gtex_by_gene <- unique(tv[, .(variant_id, Gene_name, gtex_slope_minor)], by = c("variant_id", "Gene_name"))
gene_dir <- lanes[sig_lane %in% c("DE_high", "DE_low"), .(Gene_name, deviation_sign = sign(norm_log2FC))]
panels <- orient_panels(panels, variant_alleles, gtex_by_gene, gene_dir)
# ALT is the plotted allele when both or neither of "ALT is minor" and
# "the panel plots the minor allele" hold.
panels[, alt_is_plot := alt_is_minor == plot_minor]
stopifnot(!anyNA(panels$alt_is_plot), all(panels$plot_allele == fifelse(panels$alt_is_plot, panels$ALT, panels$REF)))
panels[, set := factor(SET_LAB[as.integer(panel_group)], levels = SET_LAB)]
cat(sprintf("  panels: %s\n", paste(sprintf("%s %d", SET_LAB, tabulate(panels$set, 4)), collapse = "; ")))

# =============================================================================
# STEP 2: T21 aFC of the plotted allele, ploidy 3
# =============================================================================

cohort   <- run$cohort()
E        <- run$expression()
meta_t21 <- cohort[Karyotype == "T21"]
meta_t21[, subject_id := subject_id_from_labid(LabID)]
# Joined by Ensembl id: GTEx gene symbols differ from the roster's for some
# lncRNAs, and every roster gene carries its stable id.
roster <- fread(run$processed("eqtl_supported_genes.csv"))[, .(Gene_name, ensembl_stable)]
eg <- fread(EGENES)[gene_chr == "chr21",
                    .(ensembl_stable = sub("\\..*$", "", gene_id), lead_variant = variant_id, lead_slope = slope,
                      lead_slope_se = slope_se, lead_afc = afc, lead_afc_se = afc_se, gtex_qval = qval)]
eg <- merge(unique(roster, by = "Gene_name"), eg, by = "ensembl_stable")[, ensembl_stable := NULL]
own <- merge(unique(panels[set != "Negative control", .(Gene_name, set)]), eg, by = "Gene_name", all.x = TRUE)

geno <- fread(run$processed("genotypes_filtered.csv"),
              select = c("variant_id", "subject_id", "karyotype", "alt_dosage"))
geno_t21 <- geno[karyotype == "T21" & !is.na(alt_dosage) & subject_id %in% meta_t21$subject_id &
                   variant_id %in% c(panels$variant_id, own$lead_variant),
                 .(alt_dosage = mean(alt_dosage)), by = .(variant_id, subject_id)]
subj <- sort(unique(geno_t21$subject_id))
genes_needed <- unique(panels$Gene_name)
expr_tbl <- data.table(Gene_name = rep(genes_needed, each = length(subj)),
                       subject_id = rep(subj, times = length(genes_needed)))
lab_for_subj <- setNames(meta_t21$LabID, meta_t21$subject_id)
expr_tbl[, expr := E[cbind(match(Gene_name, rownames(E)), match(lab_for_subj[subject_id], colnames(E)))]]
cat(sprintf("  T21 subjects: %d; genes: %d\n", length(subj), length(genes_needed)))

#' aFC fits for a table of (Gene_name, variant_id, flip) rows; `flip` TRUE
#' counts the REF allele instead of ALT.
fit_pairs <- function(pairs) {
  long <- merge(pairs, geno_t21, by = "variant_id", allow.cartesian = TRUE)
  long <- merge(long, expr_tbl, by = c("Gene_name", "subject_id"))[!is.na(expr)]
  long[, dosage := fifelse(flip, T21_CHR21_PLOIDY - alt_dosage, alt_dosage)]
  long[, {
    f <- fit_afc(dosage, expr, T21_CHR21_PLOIDY)
    l <- fit_variants(matrix(dosage, ncol = 1), expr)
    list(t21_n = f$n, t21_afc = f$afc, t21_afc_se = f$se, t21_at_bound = f$at_bound,
         t21_slope = l$slope[1], t21_slope_se = l$se[1], t21_p = l$p[1])
  }, by = .(key, Gene_name, variant_id)]
}

panels[, key := sprintf("%s|%s", set, Gene_name)]
t21 <- fit_pairs(panels[, .(key, Gene_name, variant_id, flip = !alt_is_plot)])
res <- merge(panels[, .(key, set, Gene_name, variant_id, plot_allele, plot_allele_role, minor_allele,
                        gtex_maf, alt_is_plot, decoy_gene, q_gene_bh, detected)],
             t21[, !c("Gene_name", "variant_id")], by = "key")

# =============================================================================
# STEP 3: GTEx aFC at the same variant, on the same allele
# =============================================================================

gtex_best <- unique(tv[, .(Gene_name, variant_id, gtex_slope = slope, gtex_slope_se = slope_se)],
                    by = c("Gene_name", "variant_id"))
res <- merge(res, gtex_best, by = c("Gene_name", "variant_id"), all.x = TRUE, sort = FALSE)
res <- merge(res, eg, by = "Gene_name", all.x = TRUE, sort = FALSE)
# Decoy pairs have no GTEx record of their own (the merge above can pick up
# the decoy variant only as a cis variant of the plotted gene, which a decoy
# by construction is not).
res[set == "Negative control", c("gtex_slope", "gtex_slope_se") := NA_real_]
stopifnot(!anyNA(res[set != "Negative control", gtex_slope]))
res[, best_is_lead := variant_id == lead_variant]
res[, c("gtex_afc_alt", "gtex_afc_se") := gtex_afc_at_variant(gtex_slope, gtex_slope_se, lead_slope,
                                                              lead_slope_se, lead_afc, lead_afc_se)]
res[, gtex_afc := orient_effect(gtex_afc_alt, alt_is_plot)]
res[set == "Negative control", `:=`(gtex_afc = 0, gtex_afc_se = 0)]   # expectation, not an estimate
res[, gtex_source := fcase(set == "Negative control", "none (outside cis window); expected 0",
                           best_is_lead, "GTEx aFC at this variant (lead)",
                           default = "GTEx slope at this variant, scaled by the lead variant's aFC")]

res[, c("diff", "diff_se", "diff_lo", "diff_hi", "z", "p_diff") :=
      effect_difference(t21_afc, t21_afc_se, gtex_afc, gtex_afc_se)]
res[set != "Negative control", q_diff := p.adjust(p_diff, "BH")]
# A decoy panel shows the strongest of its decoy set's variants, so its naive
# z test against 0 is biased by that selection. The gene-level permutation q
# (script 03), which repeats the selection under the null, is the valid test
# of a decoy effect against 0. The naive interval is kept and plotted: it shows
# how large an aFC best-variant selection produces with no cis effect at all.
res[set == "Negative control", q_diff := q_gene_bh]
res[, `:=`(t21_lo = t21_afc - Z95 * t21_afc_se, t21_hi = t21_afc + Z95 * t21_afc_se,
           gtex_lo = gtex_afc - Z95 * gtex_afc_se, gtex_hi = gtex_afc + Z95 * gtex_afc_se,
           same_sign = sign(t21_afc) == sign(gtex_afc), ratio = t21_afc / gtex_afc)]
res[set == "Negative control", c("same_sign", "ratio") := .(NA, NA_real_)]

# Deviating genes with no GTEx variant keep a row so all of them are accounted for.
untestable <- lanes[sig_lane %in% c("DE_high", "DE_low") & !Gene_name %in% res[set %in% SET_LAB[1:2], Gene_name],
                    .(Gene_name, set = factor(lane_block(sig_lane), levels = SET_LAB),
                      status = "not testable: no GTEx cis variant")]
res[, status := "tested"]
out <- rbind(res, untestable, fill = TRUE)
setorder(out, set, Gene_name)
fwrite(out[, .(set, Gene_name, status, variant_id, plot_allele, plot_allele_role, gtex_maf, decoy_gene,
               t21_n, t21_afc, t21_afc_se, t21_lo, t21_hi, t21_at_bound, t21_slope, t21_slope_se, t21_p,
               gtex_afc, gtex_afc_se, gtex_lo, gtex_hi, gtex_source, best_is_lead, lead_variant, gtex_qval,
               diff, diff_se, diff_lo, diff_hi, p_diff, q_diff, same_sign, ratio,
               q_gene_bh, detected)], run$table("t21_vs_gtex_afc.csv"))

cat("\n  aFC per copy of the boxplot allele (log2; diff = T21 - GTEx):\n")
print(res[order(set, Gene_name), .(set, Gene_name, allele = plot_allele, lead = best_is_lead,
                                   t21 = round(t21_afc, 2), t21_se = round(t21_afc_se, 2),
                                   gtex = round(gtex_afc, 2), gtex_se = round(gtex_afc_se, 2),
                                   diff = round(diff, 2), q = signif(q_diff, 2))])
cat(sprintf("  not testable (no GTEx cis variant): %s\n", paste(untestable$Gene_name, collapse = ", ")))

# =============================================================================
# STEP 4: Sensitivity at the GTEx lead variant (no slope-to-aFC conversion)
# =============================================================================

lead_pairs <- own[!is.na(lead_variant) & lead_variant %in% geno_t21$variant_id,
                  .(key = sprintf("%s|%s", set, Gene_name), Gene_name, variant_id = lead_variant, flip = FALSE)]
lead <- merge(fit_pairs(lead_pairs), own, by = "Gene_name")
lead[, c("diff", "diff_se", "diff_lo", "diff_hi", "z", "p_diff") :=
       effect_difference(t21_afc, t21_afc_se, lead_afc, lead_afc_se)]
lead[, q_diff := p.adjust(p_diff, "BH")]
setorder(lead, set, Gene_name)
fwrite(lead[, .(set, Gene_name, lead_variant, allele = "ALT", t21_n, t21_afc, t21_afc_se, t21_at_bound,
                gtex_afc = lead_afc, gtex_afc_se = lead_afc_se, gtex_qval, diff, diff_se, diff_lo,
                diff_hi, p_diff, q_diff)], run$table("t21_vs_gtex_afc_lead_variant.csv"))

# =============================================================================
# STEP 5: Agreement summary
# =============================================================================

#' Agreement of T21 with GTEx over a set of panels; the through-origin slope
#' is inverse-variance weighted on the difference SE.
agreement <- function(d, label, x = "gtex_afc") {
  d <- d[!is.na(d$t21_afc) & !is.na(d[[x]]) & !is.na(d$diff_se)]
  w <- 1 / d$diff_se^2
  data.table(group = label, n = nrow(d),
             same_sign = sum(sign(d$t21_afc) == sign(d[[x]])),
             ci_of_difference_covers_0 = sum(d$diff_lo <= 0 & d$diff_hi >= 0),
             differs_at_fdr_05 = sum(d$q_diff < ALPHA),
             pearson_r = if (nrow(d) > 2) cor(d$t21_afc, d[[x]]) else NA_real_,
             slope_t21_on_gtex = sum(w * d$t21_afc * d[[x]]) / sum(w * d[[x]]^2),
             median_ratio = median(d$t21_afc / d[[x]]),
             median_abs_diff = median(abs(d$diff)))
}
cis <- res[set != "Negative control"]
summary_tbl <- rbind(
  agreement(cis, "All cis panels, best variant"),
  agreement(cis[set %in% SET_LAB[1:2]], "Deviating genes, best variant"),
  agreement(cis[set == "Positive control"], "Positive controls, best variant"),
  agreement(lead, "All cis genes, GTEx lead variant", x = "lead_afc"),
  res[set == "Negative control",
      .(group = "Negative controls (expected 0; FDR from permutation q)", n = .N, same_sign = NA_integer_,
        ci_of_difference_covers_0 = sum(diff_lo <= 0 & diff_hi >= 0),
        differs_at_fdr_05 = sum(q_diff < ALPHA), pearson_r = NA_real_,
        slope_t21_on_gtex = NA_real_, median_ratio = NA_real_, median_abs_diff = median(abs(diff)))])
fwrite(summary_tbl, run$table("t21_vs_gtex_afc_summary.csv"))
cat("\n  Agreement:\n")
print(summary_tbl[, .(group, n, same_sign, ci_covers_0 = ci_of_difference_covers_0,
                      fdr_05 = differs_at_fdr_05, r = round(pearson_r, 3),
                      slope = round(slope_t21_on_gtex, 2), med_ratio = round(median_ratio, 2),
                      med_abs_diff = round(median_abs_diff, 2))])

# =============================================================================
# STEP 6: Figure
# =============================================================================

theme_fig <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 9.5),
        plot.subtitle = element_text(size = 7.5, colour = "grey30"), legend.position = "bottom",
        legend.text = element_text(size = 7.5), legend.margin = margin(0, 0, 0, 0))
all_row <- summary_tbl[group == "All cis panels, best variant"]

# A: T21 against GTEx, one point per cis panel; negative controls at GTEx = 0.
# Axis limits come from the cis panels; a decoy on a rare variant can carry an
# interval many times wider, which is clipped at the panel edge rather than
# allowed to set the scale (WIDE_SE marks those rows in panel B).
WIDE_SE <- 1
cis_rows <- res[set != "Negative control"]
lim <- range(c(cis_rows$t21_lo, cis_rows$t21_hi, cis_rows$gtex_lo, cis_rows$gtex_hi), na.rm = TRUE) * 1.06
pa <- ggplot(res, aes(x = gtex_afc, y = t21_afc, colour = set, shape = set)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40", linewidth = 0.4, linetype = "dashed") +
  geom_errorbar(aes(ymin = t21_lo, ymax = t21_hi), width = 0, linewidth = 0.35, alpha = 0.7) +
  geom_errorbar(data = res[set != "Negative control"], aes(xmin = gtex_lo, xmax = gtex_hi),
                orientation = "y", width = 0, linewidth = 0.35, alpha = 0.7) +
  geom_point(size = 2.2) +
  geom_text_repel(data = res[set != "Negative control"], aes(label = Gene_name), size = 2.4,
                  fontface = "italic", colour = "grey15", min.segment.length = 0.2,
                  segment.colour = "grey60", segment.size = 0.25, box.padding = 0.3,
                  max.overlaps = Inf, seed = 1) +
  scale_colour_manual(NULL, values = COL_SET, drop = FALSE) +
  scale_shape_manual(NULL, values = SHP_SET, drop = FALSE) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(tag = "A", title = "Allelic fold change: T21 vs GTEx whole blood",
       subtitle = sprintf("Per copy of each boxplot's x-axis allele; T21 fitted at ploidy 3. %d cis panels: r = %.2f, median T21/GTEx = %.2f",
                          all_row$n, all_row$pearson_r, all_row$median_ratio),
       x = expression("GTEx whole blood"~log[2]~"aFC (euploid)"),
       y = expression("T21"~log[2]~"aFC")) +
  theme_fig

# B: the difference per panel, with the allele named.
res[, row_lab := sprintf("%s  (%s)", Gene_name, plot_allele)]
res[, flag := fcase(is.na(q_diff), "", q_diff < ALPHA, "*", default = "")]
setorder(res, set, -diff)
res[, row_id := factor(key, levels = rev(key))]
narrow <- res[is.na(diff_se) | diff_se < WIDE_SE]
xlim_b <- range(c(narrow$diff_lo, narrow$diff_hi, 0), na.rm = TRUE) + c(-0.15, 0.35)
res[, wide := !is.na(diff_se) & diff_se >= WIDE_SE]
res[, diff_shown := pmin(pmax(diff, xlim_b[1]), xlim_b[2])]
pb <- ggplot(res, aes(y = row_id, x = diff, colour = set, shape = set)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4, linetype = "dashed") +
  geom_errorbar(aes(xmin = diff_lo, xmax = diff_hi), orientation = "y", width = 0, linewidth = 0.45) +
  geom_point(data = res[wide == FALSE], size = 2.2) +
  geom_text(data = res[wide == TRUE], aes(x = 0, label = sprintf("%.1f, CI %.1f to %.1f (off scale)", diff, diff_lo, diff_hi)),
            colour = "grey25", size = 2.3, vjust = -0.7) +
  geom_text(aes(x = xlim_b[2], label = flag), colour = "grey15", size = 3.5, hjust = 1, vjust = 0.75) +
  facet_grid(set ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = setNames(res$row_lab, res$key)) +
  coord_cartesian(xlim = xlim_b) +
  scale_colour_manual(values = COL_SET, guide = "none") +
  scale_shape_manual(values = SHP_SET, guide = "none") +
  labs(tag = "B", title = "Difference, T21 minus GTEx (95% CI)",
       subtitle = "Gene (plotted allele). * differs at FDR 0.05. Negative controls: vs 0, permutation q",
       x = expression("Difference in"~log[2]~"aFC"), y = NULL) +
  theme_fig + theme(axis.text.y = element_text(face = "italic", size = 7),
                    panel.grid.major.y = element_blank(),
                    strip.text.y = element_text(angle = 0, size = 7.5, hjust = 0),
                    strip.background = element_rect(fill = "grey93", colour = NA))

# C: spike-in power per deviating gene (script 15), D: base rate of the
# detection result on Expected-dosage genes (script 16). Both are read from
# those scripts' tables; when either is missing the figure keeps A and B only.
power_path <- run$table("spike_in_power.csv"); genes_path <- run$table("spike_in_power_genes.csv")
base_path  <- run$table("base_rate_summary.csv"); pos_path <- run$table("eqtl_control_positive.csv")
have_cd <- all(file.exists(c(power_path, genes_path, base_path, pos_path)))
if (!have_cd) cat("  spike-in power or base-rate tables missing (scripts 15, 16): panels C and D skipped\n")
if (have_cd) {
  AFC_GRID <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2)
  pw <- fread(power_path); pg <- fread(genes_path)
  dirs <- lanes[sig_lane %in% c("DE_high", "DE_low"), .(Gene_name, direction = fifelse(norm_log2FC > 0, "up", "down"))]
  pw <- merge(pw, dirs, by = "Gene_name"); pg <- merge(pg[!is.na(gtex_afc)], dirs, by = "Gene_name")
  pg[, marker := fifelse(detected_observed, "detected in T21", "not detected in T21")]
  pc <- ggplot(pw, aes(x = afc, y = power_p05, group = Gene_name, colour = direction)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
    geom_line(linewidth = 0.5, alpha = 0.8) +
    geom_point(data = pg, aes(x = abs(gtex_afc), y = power_at_gtex_afc_p05, shape = marker), size = 2.6, fill = "white", stroke = 1) +
    geom_text_repel(data = pg, aes(x = abs(gtex_afc), y = power_at_gtex_afc_p05, label = Gene_name), size = 2.5,
                    fontface = "italic", colour = "grey15", max.overlaps = Inf, seed = 1, box.padding = 0.3) +
    scale_colour_manual(values = COL_DIR, guide = "none") +
    scale_shape_manual(NULL, values = c("detected in T21" = 16, "not detected in T21" = 21)) +
    scale_x_continuous(breaks = AFC_GRID) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(tag = "C", title = "Spike-in power of the test, per deviating gene",
         subtitle = "Own genotypes, variant set and noise; point: power at the gene's GTEx aFC; dashed: 80%",
         x = expression("Spiked"~log[2]~"aFC (ploidy 3)"), y = "Power (gene-level p < 0.05)") +
    theme_fig
  # D: detection rates with Wilson 95% intervals
  br <- fread(base_path); pos <- fread(pos_path)
  wilson <- function(k, n, z = Z95) { p <- k / n; d <- 1 + z^2 / n
    c <- (p + z^2 / (2 * n)) / d; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d; c(c - h, c + h) }
  rates <- rbind(
    data.table(set = "Deviating genes", k = br$n_detected[1], n = br$n_tested[1]),
    data.table(set = "Matched positive controls", k = sum(pos$detected), n = sum(!is.na(pos$p_gene_perm))),
    data.table(set = "Expected-dosage genes", k = br$n_detected[2], n = br$n_tested[2]),
    data.table(set = "Expected-dosage genes,\nsame expression range", k = br$n_detected[3], n = br$n_tested[3]))
  ci <- t(mapply(wilson, rates$k, rates$n)); rates[, `:=`(rate = k / n, lo = ci[, 1], hi = ci[, 2])]
  rates[, set := factor(set, levels = rev(set))]
  rates[, lab := sprintf("%d / %d", k, n)]
  pd <- ggplot(rates, aes(y = set, x = rate)) +
    geom_vline(xintercept = rates[set == "Deviating genes", rate], linetype = "dashed", colour = "grey55", linewidth = 0.35) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.25, colour = "grey30", linewidth = 0.45) +
    geom_point(size = 3, colour = "grey15") +
    geom_text(aes(label = lab), vjust = -1.1, size = 2.6, colour = "grey25") +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(tag = "D", title = "How often the test detects a cis-eQTL",
         subtitle = sprintf("Same test, same cohort; Wilson 95%% CI. Deviating vs expected-dosage: Fisher p = %.2f",
                            br$fisher_p_vs_deviating[2]),
         x = "Genes with a detected cis-eQTL", y = NULL) +
    theme_fig + theme(panel.grid.major.y = element_blank())
}

fig <- if (have_cd) (pa + pb + plot_layout(widths = c(1.15, 1))) / (pc + pd + plot_layout(widths = c(1.3, 1))) +
                      plot_layout(heights = c(1.15, 1)) else pa + pb + plot_layout(widths = c(1.15, 1))

save_fig <- function(p, stem, width, height) {
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height, limitsize = FALSE)
  ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = 200, limitsize = FALSE)
  for (f in paste0(stem, c(".pdf", ".png"))) {
    if (!file.exists(f)) stop("failed to write ", f)
    cat(sprintf("  Saved: %s\n", f))
  }
}
cat("\n")
save_fig(fig, run$figure("t21_vs_gtex_allelic_effect"), width = 13, height = if (have_cd) 12.5 else 6.6)

writeLines(capture.output(sessionInfo()), run$table("t21_vs_gtex_afc_session_info.txt"))
cat("\nDone.\n")
