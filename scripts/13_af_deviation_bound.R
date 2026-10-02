# 13_af_deviation_bound.R
#
# Purpose: Test whether the cis-eQTLs of the deviating chr21 genes can account
#          for the deviations themselves (manuscript review, concern 1).
#
#   The within-T21 eQTL test asks which people with T21 express a gene more.
#   A deviation is a group property: the ploidy-corrected T21/D21 ratio. Under
#   an additive allelic model a cis-eQTL moves the group ratio only through an
#   allele-frequency difference between the groups (scripts/lib/af_bound.R),
#   so for every deviating gene's cis variants this script reports
#
#     1. the ALT frequency in the run's T21 subjects (subject bootstrap CI)
#        against external references (GTEx whole blood, gnomAD v4.1 global
#        and non-Finnish European) and the 14 genotyped HTP controls;
#     2. the shift of the group log2 ratio that the observed frequency
#        difference predicts, with a CI from the slope SE and the bootstrap;
#     3. the reachable bound: the largest shift the variant could produce if
#        its ALT allele were absent from, or fixed in, T21 (at the 95% CI edge
#        of the slope, so the bound is generous);
#     4. the T21 frequency each variant would need to produce the whole
#        deviation alone, and the easiest such route over every common
#        (MAF >= 0.05), within-T21-supported (p < 0.05) cis variant, set
#        against the |delta AF| the cohort shows at common variants;
#     5. the within-T21 R-squared, which is what the eQTL does explain.
#
#   The primary reference is the external one closest to the T21 cohort over
#   every tested cis variant (smallest median |delta AF|), a data-driven proxy
#   for ancestry matching; the others are kept as columns. The background
#   distribution of delta AF over all tested variants also shows any
#   systematic genotyping offset (e.g. from merging with --missing-to-ref).
#
# Inputs (per run):
#   - tables/chr21_lane_assignments.csv (deviating genes, best variant)
#   - tables/deseq2_all_genes_ploidy_normalized.csv (lfcSE)
#   - tables/t21_dosage_per_variant.csv (within-T21 fits, GTEx/gnomAD AF)
#   - processed/genotypes_filtered.csv, processed/analysis_cohort.csv
#   - data/chr21_ctrl_PASS.csv (the 14 genotyped HTP controls)
# Outputs (per run):
#   - tables/af_deviation_bound_variants.csv   every (gene, cis variant)
#   - tables/af_deviation_bound_genes.csv      per gene: best variant + max over variants
#   - tables/af_reference_background.csv       delta AF summary per reference
#   - figures/af_deviation_bound.{pdf,png}
#
# Usage: Rscript scripts/13_af_deviation_bound.R --run adjusted

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})
source("scripts/lib/alleles.R")   # T21_CHR21_PLOIDY
source("scripts/lib/af_bound.R")
source("scripts/lib/run.R"); run <- load_run()

N_BOOT        <- 2000
Z95           <- qnorm(0.975)
MIN_MAF       <- 0.05   # "common" in the reference, for the easiest-route search
ALPHA_SUPPORT <- 0.05   # within-T21 nominal p for a variant to count as an eQTL here
COL_DIR <- c(up = "#B2182B", down = "#2166AC")
REF_LABELS <- c(gtex_af = "GTEx whole blood", gnomad_af = "gnomAD v4.1 global",
                gnomad_af_nfe = "gnomAD v4.1 NFE", ctrl_af = "HTP controls (n = 14)")

cat(sprintf("=== T21-eQTL: allele-frequency bound on the deviations [run: %s] ===\n\n", run$name))

# =============================================================================
# STEP 1: Deviating genes, observed deviation and its SE
# =============================================================================

lanes <- fread(run$table("chr21_lane_assignments.csv"))
de    <- fread(run$table("deseq2_all_genes_ploidy_normalized.csv"))
dev <- lanes[sig_lane %in% c("DE_high", "DE_low") & eqtl_lane %in% c("cis_eqtl", "no_cis_eqtl"),
             .(EnsemblID, Gene_name, sig_lane, tier, eqtl_lane, q_gene_bh,
               deviation = norm_log2FC, best_variant)]
dev <- merge(dev, de[, .(EnsemblID, de_log2FC = log2FoldChange, lfcSE)], by = "EnsemblID")
stopifnot(nrow(dev) == sum(lanes$eqtl_lane %in% c("cis_eqtl", "no_cis_eqtl")),
          all(abs(dev$de_log2FC - dev$deviation) < 1e-8))
dev[, `:=`(dev_lo = deviation - Z95 * lfcSE, dev_hi = deviation + Z95 * lfcSE, de_log2FC = NULL)]
cat(sprintf("  Tested deviating genes: %d\n", nrow(dev)))

# =============================================================================
# STEP 2: Within-T21 fits and reference frequencies for every cis variant
# =============================================================================

fits <- fread(run$table("t21_dosage_per_variant.csv"))
fits <- fits[Gene_name %in% dev$Gene_name & !is.na(t21_slope),
             .(Gene_name, variant_id, POS, REF, ALT, t21_n, t21_slope, t21_se, t21_t, t21_p,
               gtex_af, gnomad_af, gnomad_af_nfe)]
cat(sprintf("  (gene, variant) fits: %d over %d variants\n", nrow(fits), uniqueN(fits$variant_id)))

# =============================================================================
# STEP 3: T21 cohort ALT frequency with a subject bootstrap
# =============================================================================

cohort <- run$cohort()
t21_subjects <- cohort[Karyotype == "T21" & (is.na(excluded_reason) | excluded_reason == ""),
                       unique(sub("[A-Z][0-9]*$", "", LabID))]
geno <- fread(run$processed("genotypes_filtered.csv"),
              select = c("variant_id", "subject_id", "karyotype", "alt_dosage"))
geno <- geno[karyotype == "T21" & subject_id %in% t21_subjects & !is.na(alt_dosage),
             .(alt_dosage = mean(alt_dosage)), by = .(variant_id, subject_id)]
D <- as.matrix(dcast(geno[variant_id %in% fits$variant_id], variant_id ~ subject_id,
                     value.var = "alt_dosage"), rownames = "variant_id")
if (anyNA(D)) stop("T21 dosage matrix has missing subject-variant cells")
cat(sprintf("  T21 subjects: %d; variants with genotypes: %d\n", ncol(D), nrow(D)))

af_boot <- bootstrap_alt_af(D, T21_CHR21_PLOIDY, B = N_BOOT, seed = 42)
t21_af <- data.table(variant_id = rownames(D),
                     t21_af    = rowMeans(D) / T21_CHR21_PLOIDY,
                     t21_af_lo = apply(af_boot, 1, quantile, 0.025),
                     t21_af_hi = apply(af_boot, 1, quantile, 0.975))

# =============================================================================
# STEP 4: The 14 genotyped HTP controls (small; a within-study check only)
# =============================================================================

ctrl_hdr  <- names(fread("data/chr21_ctrl_PASS.csv", nrows = 0))
ctrl_ids  <- setdiff(ctrl_hdr, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"))
ctrl <- fread("data/chr21_ctrl_PASS.csv", select = c("POS", "REF", "ALT", ctrl_ids))
ctrl <- ctrl[POS %in% fits$POS]
ctrl_gt <- as.matrix(ctrl[, ..ctrl_ids])
ctrl_gt <- sub(":.*$", "", ctrl_gt)                              # GT field only
ctrl_alt <- matrix(nchar(gsub("[^1-9]", "", ctrl_gt)), nrow(ctrl_gt))   # ALT copies
ctrl[, ctrl_af := rowSums(ctrl_alt) / (2 * length(ctrl_ids))]
ctrl[, variant_id := sprintf("chr21_%d_%s_%s_b38", POS, REF, ALT)]
# --missing-to-ref merging: a variant absent from the control call set is REF
# in all 14, matching how the T21 set treats absent calls.
ctrl_af <- ctrl[!duplicated(variant_id), .(variant_id, ctrl_af)]

# =============================================================================
# STEP 5: Background delta AF over every tested variant; primary reference
# =============================================================================

vars <- unique(fits[, .(variant_id, gtex_af, gnomad_af, gnomad_af_nfe)])
vars <- merge(vars, t21_af, by = "variant_id")
vars <- merge(vars, ctrl_af, by = "variant_id", all.x = TRUE)
vars[is.na(ctrl_af), ctrl_af := 0]

ref_cols <- names(REF_LABELS)
background <- rbindlist(lapply(ref_cols, function(rc) {
  d <- vars$t21_af - vars[[rc]]
  data.table(reference = rc, label = REF_LABELS[[rc]], n_variants = sum(!is.na(d)),
             median_delta = median(d, na.rm = TRUE),
             median_abs_delta = median(abs(d), na.rm = TRUE),
             q95_abs_delta = quantile(abs(d), 0.95, na.rm = TRUE),
             cor = cor(vars$t21_af, vars[[rc]], use = "complete.obs"))
}))
PRIMARY <- background[reference != "ctrl_af"][which.min(median_abs_delta), reference]
cat("\n  Background delta AF (T21 minus reference) over tested cis variants:\n")
print(background[, .(label, n_variants, median_delta = signif(median_delta, 3),
                     median_abs_delta = signif(median_abs_delta, 3),
                     q95_abs_delta = signif(q95_abs_delta, 3), cor = round(cor, 4))])
cat(sprintf("  Primary reference (closest to the T21 cohort): %s\n", REF_LABELS[[PRIMARY]]))
fwrite(background[, primary := reference == PRIMARY], run$table("af_reference_background.csv"))

# =============================================================================
# STEP 6: Predicted shift, reachable bound and within-T21 R2 per (gene, variant)
# =============================================================================

vt <- merge(fits[, .(Gene_name, variant_id, t21_n, t21_slope, t21_se, t21_t, t21_p)],
            vars, by = "variant_id")
vt <- merge(vt, dev[, .(Gene_name, deviation, best_variant)], by = "Gene_name")
vt[, p_ref := get(PRIMARY)]
vt <- vt[!is.na(p_ref)]
vt[, `:=`(delta_af = t21_af - p_ref,
          r = allelic_r_from_slope(t21_slope, T21_CHR21_PLOIDY),
          r2_within_t21 = r2_from_t(t21_t, t21_n))]
vt[, pred_shift := predicted_group_shift(r, t21_af, p_ref)]

# Uncertainty of the predicted shift: slope drawn from its sampling
# distribution, T21 frequency from the subject bootstrap (reference fixed).
set.seed(7)
slope_draw <- matrix(rnorm(nrow(vt) * N_BOOT, vt$t21_slope, vt$t21_se), nrow(vt))
r_draw     <- allelic_r_from_slope(slope_draw, T21_CHR21_PLOIDY)
af_draw    <- af_boot[match(vt$variant_id, rownames(D)), , drop = FALSE]
shift_draw <- log2((1 + r_draw * af_draw) / (1 + r_draw * vt$p_ref))
vt[, `:=`(pred_lo = apply(shift_draw, 1, quantile, 0.025),
          pred_hi = apply(shift_draw, 1, quantile, 0.975))]

# Reachable bound in the deviation's direction, generous: the larger of the
# bounds at the two 95% CI edges of the slope.
edge_lo <- shift_bounds(allelic_r_from_slope(vt$t21_slope - Z95 * vt$t21_se), vt$p_ref)
edge_hi <- shift_bounds(allelic_r_from_slope(vt$t21_slope + Z95 * vt$t21_se), vt$p_ref)
point   <- shift_bounds(vt$r, vt$p_ref)
vt[, `:=`(bound_low = pmin(edge_lo$bound_low, edge_hi$bound_low),
          bound_high = pmax(edge_lo$bound_high, edge_hi$bound_high),
          bound_point_toward = bound_toward(point$bound_low, point$bound_high, deviation))]
vt[, bound_toward_dev := bound_toward(bound_low, bound_high, deviation)]
# The T21 frequency each variant would need to produce the whole deviation on
# its own (point slope), as a difference from the reference; NA = out of reach.
vt[, `:=`(req_af = required_case_af(r, p_ref, deviation),
          frac_pred = pred_shift / deviation, frac_bound = bound_toward_dev / deviation,
          is_best = variant_id == best_variant,
          common = pmin(p_ref, 1 - p_ref) >= MIN_MAF,
          supported = t21_p < ALPHA_SUPPORT)]
vt[, req_delta_af := req_af - p_ref]
setorder(vt, Gene_name, t21_p)
fwrite(vt[, !c("best_variant")], run$table("af_deviation_bound_variants.csv"))

# Background for a required difference: the |delta AF| the cohort actually
# shows at common variants, whatever their eQTL status. The figure shades the
# 99th percentile rather than the 99.9th: the extreme tail is a few indel
# records whose common ALT allele went uncalled in T21 and was merged as REF
# (--missing-to-ref), e.g. the chr21:14.23 Mb cluster near RBM11 in baseline.
bg_common <- abs(vars[pmin(get(PRIMARY), 1 - get(PRIMARY)) >= MIN_MAF, t21_af - get(PRIMARY)])
BG_Q <- quantile(bg_common, c(0.99, 0.999), na.rm = TRUE)
cat(sprintf("\n  Background |delta AF| at common variants (MAF >= %.2f, n = %d): 99th %.3f, 99.9th %.3f, max %.3f\n",
            MIN_MAF, sum(!is.na(bg_common)), BG_Q[1], BG_Q[2], max(bg_common, na.rm = TRUE)))

# Per gene: the best variant of the gene-level test, plus the easiest route
# any common, within-T21-supported cis variant offers (smallest required
# |delta AF|; out-of-reach variants count as impossible).
best <- vt[is_best == TRUE, .(Gene_name, variant_id, t21_slope, t21_se, r2_within_t21,
                              t21_af, t21_af_lo, t21_af_hi, p_ref, delta_af, ctrl_af,
                              gtex_af, gnomad_af, gnomad_af_nfe,
                              pred_shift, pred_lo, pred_hi, bound_low, bound_high,
                              bound_toward_dev, frac_pred, frac_bound,
                              req_af, req_delta_af)]
easiest <- vt[common & supported, {
  i <- if (all(is.na(req_delta_af))) NA_integer_ else which.min(abs(req_delta_af))
  .(n_variants_common_supported = .N,
    easiest_variant = if (is.na(i)) NA_character_ else variant_id[i],
    easiest_req_delta_af = if (is.na(i)) NA_real_ else req_delta_af[i],
    easiest_obs_delta_af = if (is.na(i)) NA_real_ else delta_af[i],
    n_reachable = sum(!is.na(req_delta_af)))
}, by = Gene_name]
genes <- merge(dev, best, by = "Gene_name", all.x = TRUE)
genes <- merge(genes, easiest, by = "Gene_name", all.x = TRUE)
genes[, `:=`(reference = REF_LABELS[[PRIMARY]],
             easiest_req_pctl = vapply(abs(easiest_req_delta_af), function(x)
               if (is.na(x)) NA_real_ else mean(bg_common < x, na.rm = TRUE), numeric(1)))]
# Deviation-aligned copies: every shift signed so that positive points the way
# the gene deviates (up for DE_high, down for DE_low). This removes the
# up/down mirror between the two groups and leaves one question per row: how
# far toward the deviation can the variant move the group ratio?
genes[, dev_sign := sign(deviation)]
genes[, `:=`(toward_obs = abs(deviation),
             toward_obs_lo = pmin(dev_sign * dev_lo, dev_sign * dev_hi),
             toward_obs_hi = pmax(dev_sign * dev_lo, dev_sign * dev_hi),
             toward_pred = dev_sign * pred_shift,
             toward_pred_lo = pmin(dev_sign * pred_lo, dev_sign * pred_hi),
             toward_pred_hi = pmax(dev_sign * pred_lo, dev_sign * pred_hi),
             toward_reach = pmax(dev_sign * bound_low, dev_sign * bound_high),
             away_reach = pmin(dev_sign * bound_low, dev_sign * bound_high))]
setorder(genes, -deviation)
fwrite(genes, run$table("af_deviation_bound_genes.csv"))

cat("\n  Per gene (log2 units; AF differences are T21 minus reference):\n")
print(genes[, .(Gene_name, eqtl_lane, deviation = round(deviation, 2),
                obs_dAF = round(delta_af, 3), pred = round(pred_shift, 3),
                pred_ci = sprintf("[%.3f, %.3f]", pred_lo, pred_hi),
                bound = round(bound_toward_dev, 2),
                req_dAF_best = round(req_delta_af, 2),
                req_dAF_easiest = round(easiest_req_delta_af, 2),
                n_reach = sprintf("%d/%d", n_reachable, n_variants_common_supported),
                r2 = round(r2_within_t21, 3))])

# =============================================================================
# STEP 7: Figure
# =============================================================================

genes[, direction := fifelse(deviation > 0, "up", "down")]
genes[, gene_lab := sprintf("%s  (%s)", Gene_name,
                            fifelse(eqtl_lane == "cis_eqtl", "cis-eQTL", "not detected"))]
genes[, gene_lab := factor(gene_lab, levels = rev(gene_lab))]

# A: T21 cohort frequency against the primary reference, every tested variant.
bg <- background[reference == PRIMARY]
pa <- ggplot(vars, aes(x = get(PRIMARY), y = t21_af)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_point(size = 0.5, alpha = 0.25, colour = "grey45", stroke = 0) +
  geom_point(data = genes, aes(x = p_ref, y = t21_af, colour = direction),
             size = 2.2, shape = 21, fill = "white", stroke = 0.9) +
  geom_text_repel(data = genes, aes(x = p_ref, y = t21_af, label = Gene_name, colour = direction),
                  size = 2.6, fontface = "italic", min.segment.length = 0, seed = 1,
                  box.padding = 0.35, max.overlaps = Inf, show.legend = FALSE) +
  scale_colour_manual(values = COL_DIR, guide = "none") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(tag = "A", title = "ALT allele frequency, T21 cohort vs reference",
       subtitle = sprintf("%s tested cis variants; median T21 - ref = %+.4f, median absolute diff = %.3f",
                          format(bg$n_variants, big.mark = ","), bg$median_delta, bg$median_abs_delta),
       x = sprintf("%s ALT frequency", REF_LABELS[[PRIMARY]]),
       y = sprintf("T21 cohort ALT frequency (n = %d)", ncol(D))) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 9.5),
        plot.subtitle = element_text(size = 7.5, colour = "grey30"))

# B: observed deviation vs what the best variant could do, on the
# deviation-aligned axis (positive = toward the gene's deviation, for up and
# down genes alike). Group ratios use 3 chr21 copies in T21 and 2 in D21.
LAB_B <- c(obs = "Observed deviation (95% CI)",
           pred = "Shift predicted by the observed AF difference (95% CI)",
           rng = "Reachable if T21 AF were 0 or 1")
theme_rows <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 7.5),
        legend.spacing.y = unit(0, "pt"), legend.margin = margin(0, 0, 0, 0),
        plot.title = element_text(face = "bold", size = 9.5),
        plot.subtitle = element_text(size = 7.5, colour = "grey30"))
pb <- ggplot(genes, aes(y = gene_lab)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_linerange(aes(xmin = away_reach, xmax = toward_reach, linewidth = LAB_B[["rng"]]),
                 colour = "grey82") +
  geom_errorbar(aes(xmin = toward_pred_lo, xmax = toward_pred_hi), orientation = "y", width = 0.25,
                colour = "black", linewidth = 0.4) +
  geom_point(aes(x = toward_pred, shape = LAB_B[["pred"]]), size = 2.2, fill = "black") +
  geom_errorbar(aes(xmin = toward_obs_lo, xmax = toward_obs_hi, colour = direction), orientation = "y",
                width = 0.25, linewidth = 0.5) +
  geom_point(aes(x = toward_obs, colour = direction, shape = LAB_B[["obs"]]), size = 2.6) +
  scale_colour_manual(values = COL_DIR, labels = c(up = "Higher than expected", down = "Lower than expected"),
                      name = NULL) +
  scale_shape_manual(NULL, values = setNames(c(16, 23), unname(LAB_B[c("obs", "pred")])),
                     breaks = unname(LAB_B[c("obs", "pred")])) +
  scale_linewidth_manual(NULL, values = setNames(4, unname(LAB_B[["rng"]]))) +
  guides(shape = guide_legend(order = 1, ncol = 1,
                              override.aes = list(colour = c("grey20", "black"), fill = c(NA, "black"))),
         linewidth = guide_legend(order = 2), colour = guide_legend(order = 3, override.aes = list(shape = 16))) +
  labs(tag = "B", title = "Deviation vs the best variant's reach",
       subtitle = "Aligned to each gene's deviation; T21 = 3 copies, D21 = 2; reach at the 95% CI edge of the slope",
       x = expression(symbol("\254")~"away"~~~~"Shift of"~log[2]~"(T21/D21) from 1.5"~~~~"toward the deviation"~symbol("\256")),
       y = NULL) +
  theme_rows + theme(axis.text.y = element_text(face = "italic"), plot.title = element_text(face = "bold", size = 9.5,
                                                                                          margin = margin(b = 2)))

# C: the T21 frequency change each gene's deviation would need.
LAB_C <- c(best = "Best variant",
           easy = "Easiest common variant with within-T21 p < 0.05")
req <- melt(genes[, .(gene_lab, best = abs(req_delta_af), easy = abs(easiest_req_delta_af))],
            id.vars = "gene_lab", variable.name = "which", value.name = "req")
req[, which := LAB_C[as.character(which)]]
# Out-of-reach genes (no frequency in [0, 1] produces the deviation) sit in
# their own column to the right of 1.
X_OOR <- 1.2
req[, `:=`(reachable = !is.na(req), x = fifelse(is.na(req), X_OOR, req))]
req[, y_nudge := fifelse(which == LAB_C[["best"]], 0.14, -0.14)]
pc <- ggplot(req, aes(y = gene_lab, x = x)) +
  annotate("rect", xmin = 0, xmax = BG_Q[[1]], ymin = -Inf, ymax = Inf, fill = "grey85") +
  geom_vline(xintercept = 1.1, colour = "grey70", linewidth = 0.3) +
  geom_point(aes(shape = which), size = 2.3, colour = "grey15",
             position = position_nudge(y = req$y_nudge)) +
  scale_shape_manual(NULL, values = setNames(c(16, 1), unname(LAB_C))) +
  scale_x_continuous(limits = c(0, 1.3), breaks = c(0, 0.25, 0.5, 0.75, 1, X_OOR),
                     labels = c("0", "0.25", "0.5", "0.75", "1", "out of\nreach")) +
  guides(shape = guide_legend(ncol = 1)) +
  labs(tag = "C", title = "AF difference needed to explain it",
       subtitle = sprintf("Shaded: 99%% of observed absolute T21 - ref differences (<= %.3f)", BG_Q[[1]]),
       x = "Absolute T21 - reference ALT frequency difference needed", y = NULL) +
  theme_rows + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# D: what the eQTL does explain: between-person variance within T21.
pd <- ggplot(genes, aes(y = gene_lab, x = 100 * r2_within_t21, fill = direction)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * r2_within_t21)), hjust = -0.15, size = 2.5, colour = "grey20") +
  scale_fill_manual(values = COL_DIR, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(tag = "D", title = "Within-T21 variance explained",
       subtitle = "Best variant, R-squared", x = "% of within-T21 variance", y = NULL) +
  theme_rows + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

fig <- pa + (pb + pc + pd + plot_layout(widths = c(2.2, 1.6, 1))) + plot_layout(widths = c(1, 2.6))

save_fig <- function(p, stem, width, height) {
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height, limitsize = FALSE)
  ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = 200, limitsize = FALSE)
  for (f in paste0(stem, c(".pdf", ".png"))) {
    if (!file.exists(f)) stop("failed to write ", f)
    cat(sprintf("  Saved: %s\n", f))
  }
}
cat("\n")
save_fig(fig, run$figure("af_deviation_bound"), width = 16, height = 5.8)

writeLines(capture.output(sessionInfo()), run$table("af_deviation_bound_session_info.txt"))
cat("\nDone.\n")
