# 15_spike_in_power.R
#
# Purpose: Power of the gene-level cis-eQTL test for each tested deviating
#          gene, by spike-in (manuscript review, concern 5). "Tested, not
#          detected" only means something if the test would have found the
#          gene's GTEx eQTL had it acted in T21. For each gene: the gene's own
#          T21 genotypes at its own cis variant set, its own expression noise
#          (the observed expression, permuted across subjects so no genotype
#          association remains), and an allelic effect of a given aFC added at
#          a causal variant (the GTEx lead variant when it is in the set,
#          otherwise the set's smallest-GTEx-p variant) under the ploidy-3
#          allelic model. The best-variant search is then scored against the
#          gene's own permutation null, exactly as in script 03.
#
#   Power is reported on a grid of aFC values and at the gene's GTEx aFC,
#   with the aFC needed for 80% power.
#
# Inputs (per run):
#   - tables/eqtl_gene_level_perm.csv, t21_dosage_per_variant.csv
#   - processed/genotypes_filtered.csv, expression_adjusted.csv,
#     analysis_cohort.csv, eqtl_supported_genes.csv
#   - data/Whole_Blood.v10.eGenes.txt.gz
# Outputs (per run):
#   - tables/spike_in_power.csv        gene x aFC power
#   - tables/spike_in_power_genes.csv  per gene: GTEx aFC, power there, aFC for 80%
#   - figures/spike_in_power.{pdf,png}
#
# Usage: Rscript scripts/15_spike_in_power.R --run adjusted

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})
source("scripts/lib/cohort.R")
source("scripts/lib/eqtl_fit.R")
source("scripts/lib/alleles.R")
source("scripts/lib/allelic_fc.R")
source("scripts/lib/run.R"); run <- load_run()

N_PERM  <- as.integer(Sys.getenv("T21_N_PERM", "1000"))
N_SIM   <- as.integer(Sys.getenv("T21_N_SIM", "200"))
AFC_GRID <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2)
ALPHA_GENE <- c(0.05, 0.01)
COL_DIR <- c(up = "#B2182B", down = "#2166AC")

cat(sprintf("=== T21-eQTL: spike-in power [run: %s] ===\n\n", run$name))

#' Smallest p across the columns of G for every column of Y (subjects x sims):
#' the vectorised fit of perm_min_p, applied to arbitrary expression columns.
min_p_columns <- function(G, Y) {
  n <- nrow(G)
  Gc <- center_cols(G); Sgg <- colSums(Gc^2); keep <- Sgg > 0
  Gc <- Gc[, keep, drop = FALSE]; Sgg <- Sgg[keep]
  Yc <- sweep(Y, 2, colMeans(Y), "-")
  Sge <- crossprod(Gc, Yc)
  See <- matrix(colSums(Yc^2), nrow = nrow(Sge), ncol = ncol(Y), byrow = TRUE)
  slope <- Sge / Sgg
  se <- sqrt(((See - Sge^2 / Sgg) / (n - 2)) / Sgg)
  P <- 2 * pt(-abs(slope / se), df = n - 2)
  P[!is.finite(P)] <- NA_real_
  apply(P, 2, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
}

perm  <- fread(run$table("eqtl_gene_level_perm.csv"))[!is.na(p_gene_perm)]
fits  <- fread(run$table("t21_dosage_per_variant.csv"))[Gene_name %in% perm$Gene_name]
roster <- fread(run$processed("eqtl_supported_genes.csv"))[, .(Gene_name, ensembl_stable)]
eg <- fread("data/Whole_Blood.v10.eGenes.txt.gz")[gene_chr == "chr21"]
eg[, ensembl_stable := sub("\\..*$", "", gene_id)]
eg <- merge(unique(roster, by = "Gene_name"), eg[, .(ensembl_stable, lead_variant = variant_id, gtex_afc = afc)],
            by = "ensembl_stable")

cohort <- run$cohort(); t21 <- cohort[Karyotype == "T21"]; t21[, subject_id := subject_id_from_labid(LabID)]
geno <- fread(run$processed("genotypes_filtered.csv"), select = c("variant_id", "subject_id", "karyotype", "alt_dosage"))
geno <- geno[karyotype == "T21" & subject_id %in% t21$subject_id & !is.na(alt_dosage) & variant_id %in% fits$variant_id,
             .(alt_dosage = mean(alt_dosage)), by = .(variant_id, subject_id)]
gwide <- dcast(geno, subject_id ~ variant_id, value.var = "alt_dosage")
subj <- gwide$subject_id; G_all <- as.matrix(gwide[, -1])
E <- run$expression()

power_rows <- rbindlist(lapply(perm$Gene_name, function(g) {
  vars <- intersect(fits[Gene_name == g, variant_id], colnames(G_all))
  G <- G_all[, vars, drop = FALSE]
  y <- expr_from_matrix(g, subj, E, t21)
  ok <- !is.na(y); G <- G[ok, , drop = FALSE]; y <- y[ok]; n <- length(y)
  lead <- eg[Gene_name == g, lead_variant]
  causal <- if (length(lead) && lead %in% vars) lead else fits[Gene_name == g][which.min(gtex_pval), variant_id]
  # Null: expression permuted (no association left), scored by the min-p null
  # of this gene's own variant set, as script 03 does.
  set.seed(7)
  y_null <- y[sample.int(n)]
  null_min_p <- perm_min_p(G, y_null, n_perm = N_PERM, seed = 11)
  # Simulations: fresh permutation of expression per simulation, plus the
  # allelic effect at the causal variant at each aFC of the grid.
  set.seed(13)
  Y0 <- matrix(y[as.vector(replicate(N_SIM, sample.int(n)))], nrow = n)
  g_causal <- G[, causal]
  rbindlist(lapply(AFC_GRID, function(a) {
    Y <- Y0 + afc_curve(g_causal, a, T21_CHR21_PLOIDY)
    mp <- min_p_columns(G, Y)
    pp <- vapply(mp, gene_level_p, numeric(1), min_p_perm = null_min_p)
    data.table(Gene_name = g, n_subjects = n, n_variants = length(vars), causal_variant = causal,
               causal_is_lead = identical(causal, lead), afc = a,
               power_p05 = mean(pp < ALPHA_GENE[1]), power_p01 = mean(pp < ALPHA_GENE[2]))
  }))
}))
fwrite(power_rows, run$table("spike_in_power.csv"))

# Per gene: power at the GTEx aFC (interpolated on the grid) and the aFC for 80%.
interp <- function(x, y, at) approx(x, y, xout = at, rule = 2)$y
inv    <- function(x, y, target) if (max(y) < target) NA_real_ else approx(y, x, xout = target, ties = "ordered")$y
genes_tbl <- merge(power_rows[, .(afc_80_power_p05 = inv(afc, power_p05, 0.8),
                                  causal_variant = causal_variant[1], causal_is_lead = causal_is_lead[1],
                                  n_variants = n_variants[1]), by = Gene_name],
                   eg[, .(Gene_name, gtex_afc)], by = "Gene_name", all.x = TRUE)
genes_tbl <- merge(genes_tbl, perm[, .(Gene_name, q_observed = q_gene_bh, detected_observed = cis_eqtl_detected)], by = "Gene_name")
genes_tbl[, power_at_gtex_afc_p05 := mapply(function(g, a) if (is.na(a)) NA_real_ else
  interp(power_rows[Gene_name == g, afc], power_rows[Gene_name == g, power_p05], abs(a)), Gene_name, gtex_afc)]
genes_tbl[, power_at_gtex_afc_p01 := mapply(function(g, a) if (is.na(a)) NA_real_ else
  interp(power_rows[Gene_name == g, afc], power_rows[Gene_name == g, power_p01], abs(a)), Gene_name, gtex_afc)]
setorder(genes_tbl, -power_at_gtex_afc_p05)
fwrite(genes_tbl, run$table("spike_in_power_genes.csv"))
cat("  Power at the gene's GTEx aFC (gene-level p < 0.05 / < 0.01), and aFC for 80% power:\n")
print(genes_tbl[, .(Gene_name, detected_observed, gtex_afc = round(abs(gtex_afc), 2),
                    power_p05 = round(power_at_gtex_afc_p05, 2), power_p01 = round(power_at_gtex_afc_p01, 2),
                    afc_for_80pct = round(afc_80_power_p05, 2), n_variants, causal_is_lead)])

# Figure: power curves, one line per gene, GTEx aFC marked.
lanes <- fread(run$table("chr21_lane_assignments.csv"))[, .(Gene_name, direction = fifelse(norm_log2FC > 0, "up", "down"))]
power_rows <- merge(power_rows, lanes, by = "Gene_name")
marks <- merge(genes_tbl[!is.na(gtex_afc)], lanes, by = "Gene_name")
p <- ggplot(power_rows, aes(x = afc, y = power_p05, group = Gene_name, colour = direction)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  geom_point(data = marks, aes(x = abs(gtex_afc), y = power_at_gtex_afc_p05), size = 2.4) +
  geom_text_repel(data = marks, aes(x = abs(gtex_afc), y = power_at_gtex_afc_p05, label = Gene_name),
                  size = 2.6, fontface = "italic", show.legend = FALSE, max.overlaps = Inf, seed = 1) +
  scale_colour_manual(values = COL_DIR, labels = c(up = "Higher than expected", down = "Lower than expected"), name = NULL) +
  scale_x_continuous(breaks = AFC_GRID) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Spike-in power of the gene-level cis-eQTL test",
       subtitle = sprintf("Each gene's own T21 genotypes, cis variant set and expression noise (n = %d subjects, %d simulations per point); point = the gene's GTEx aFC",
                          max(power_rows$n_subjects), N_SIM),
       x = expression("Spiked allelic fold change ("*log[2]*" aFC, ploidy 3)"),
       y = "Power (gene-level permutation p < 0.05)") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 9.5), plot.subtitle = element_text(size = 7.5, colour = "grey30"))
stem <- run$figure("spike_in_power")
ggsave(paste0(stem, ".pdf"), p, width = 7, height = 5); ggsave(paste0(stem, ".png"), p, width = 7, height = 5, dpi = 200)
figs <- paste0(stem, c(".pdf", ".png"))
if (!all(file.exists(figs))) stop("failed to write: ", paste(figs[!file.exists(figs)], collapse = ", "))
cat(sprintf("  Saved: %s.{pdf,png}\n", stem))
cat("\nDone.\n")
