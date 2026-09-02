# Decisions, legacy notes, and gotchas

The README documents the current pipeline only. This file records how it
got there: design decisions and retired approaches, legacy data sources
and terminology, the archived script catalog, and practical gotchas.

## Decision log

### Classification threshold: Hunter's rule; cohort-SD filter retired

The classification rule matches Hunter et al.'s own criterion directly:
`norm_padj < ALPHA_DE` (0.01) AND `abs(norm_log2FC) >= DEVIATION_LFC`
(`log2(1.5)`), applied on the ploidy-corrected scale. An earlier pipeline
instead gated deviations on a cohort-derived noise threshold (1 SD of the
non-chr21 cohort noise; the constant was then named `MAGNITUDE_THRESHOLD`).
That filter was retired: ploidy normalization only acts on chr21 genes, so
a non-chr21 cohort-noise SD is not a valid reference for the chr21 null,
and retiring the filter removed that mismatch. `DEVIATION_LFC` is the one
dial in this analysis.

### chr21-internal outlier test: annotation only

A chr21-internal FDR-controlled robust outlier test
(`scripts/lib/chr21_threshold.R`: median/MAD null estimated from
expressed, non-repeat chr21 genes, BH-FDR on the robust z-score at
`OUTLIER_FDR = 0.10`) is computed alongside classification and written as
annotation columns `dev_z` / `q_outlier` (plus `chr21_k_sensitivity.csv`).
It does not gate classification - only Hunter's rule does.

### eQTL classification: gene-level permutation test replaced the any-variant rule

The locus-level "any cis variant matches direction and reproduces in T21"
rule (script 03's `strongest_supp_variant` logic) is retained as context
columns (`n_cis_total`, `n_dir_match`, `n_supp_with_repro`) but is **not**
the classification rule: `results/tables/eqtl_negative_controls.csv` shows
it returns 100% "explained" even when the observed deviation direction is
artificially flipped, so it does not discriminate real signal from chance.
The classification rule is the gene-level permutation test
(`scripts/lib/eqtl_fit.R`, run in script 03): for each deviating gene, the
best-variant test statistic is compared against its null distribution
under permutation of genotype-to-expression assignment, giving
`p_gene_perm`; BH-adjusted across deviating genes to `q_gene_bh`, with
`eqtl_lane = cis_eqtl` when `q_gene_bh < FDR_GENE` (0.05).

### `eqtl_lane` is a detection result, not an "explained by eQTL" claim

`cis_eqtl` means a cis-eQTL is detectable for the gene at
`q_gene_bh < 0.05` in the within-T21 data - not that the eQTL
quantitatively accounts for the observed deviation.

### Composition-control null: correlation-matched, not independent

The composition-control null must be matched, not independent.
`partner_null(L_ctrl, lfc, n_partners = 20, n_draw = 300, seed = 1)` draws
`n_draw` random seed genes from the same non-chr21 expressed pool and, for
each seed, takes the median log2FC of *its own* top-`n_partners`
correlated genes - a null draw is a co-expression module built exactly
like the observed one. The seed is excluded from its own partner set. An
earlier version drew independent random 20-gene sets, which is
anti-conservative: a module moves together, so its median log2FC is
several times more variable (measured here: SD 0.259 matched vs 0.058
independent), and an independent null called almost any partner shift
significant. `partner_p()` reports the one-sided empirical p as
`(1 + k) / (n_draw + 1)`, so it is never exactly 0.

## Legacy GTEx source

`data/Whole_Blood.v10.eQTLs.signif_pairs.parquet` (per-gene FDR-passing
pairs only). Older runs of script 02 used this; the current pipeline reads
the chr21 allpairs extract and applies `pval_nominal <= 1e-4`
(`GTEX_PVAL_KEEP`, matching the effective signif_pairs cutoff) to keep the
variant universe manageable.

## Legacy terminology: `Sig_high_FC` vs `DE_high`

An early version of the pipeline used `Sig_high_FC` in the
`eqtl_supported_genes.csv` `gene_set` column (written by script 02). The
current lane terminology (script 04 onward) is `DE_high`. Both refer to
the same set: genes with `norm_padj < 0.01` AND `norm_log2FC >=
log2(1.5)`. The cut is on the PLOIDY-CORRECTED log2FC and split by its
sign, not on `raw_FC` - the raw chr21 fold change centres on 1.5 by ploidy
alone, so a raw-FC cut would select on the trisomy itself. Figure labels
that said "raw FC" were wrong and were corrected.

## Archived / supplementary scripts (`scripts/archive/`, do not edit)

Original paper-style Panel D pipeline:
- `02_categorize_genes` - paper-style gene categorization
- `03_volcano_plot` - diagnostic volcano
- `04_alluvial_plot` - original Panel D Sankey
- `05_eqtl_analysis` - original eQTL cross-reference (template-based)
- `06_alluvial_with_eqtl` - enhanced Panel D with eQTL terminals
- `07_sankeymatic_export` - SankeyMATIC export from old categorization
- `08_expected_dosage_eqtl` - eQTL analysis for >=1.5 FC genes

Supplementary outputs from the cohort-scale chain (run on demand):
- `10_eqtl_genotype_concordance` - per-variant T21 vs Control dosage
  means; directional concordance (`results/tables/eqtl_genotype_concordance_*.csv`)
- `14_dosage_lane_boxplots` - per-quadrant focused boxplot PDFs
  (`chr21_dosage_de_*.pdf`); consumes `t21_representative_variants.csv`,
  the legacy per-gene strongest-supportive-variant table still written by
  script 03
- `16_chr21_de_forest_plot` - per-DE-gene log2FC + 95% CI
  (`chr21_de_forest_plot.{pdf,png}`)
- `17_chr21_quadrant_plot` - within-T21 slope vs deviation scatter
  (`chr21_quadrant_plot.{pdf,png}`)

Exploratory / one-off scripts, kept but not expected to run cleanly
against the current state:
- `diagnostic_check`, `investigate_pc2`, `pca_chr21_only`,
  `process_blacklist`, `run_all.sh` (the old shell driver)

## Gotchas

- **`passes_magnitude_filter` is not "expected dosage"**: it is FALSE for
  repeat-flagged and low-expression genes too, because they are never
  eligible for the cut. Split figures on `sig_lane`, never on that flag -
  doing the latter drew 41 unassessable genes inside the Expected-dosage
  stratum in script 05 until it was fixed.
- **Self-loops in SankeyMATIC**: `chr21_lane_sankeymatic_input.txt` skips
  level-to-level passes where the source and target name are identical
  (e.g., Expected dosage genes terminate at level 2 and would otherwise
  self-loop at levels 3 and 4). The flow file
  `chr21_lane_alluvial_flow.csv` keeps them.
- **Two T21 expression-only subjects**: 304 T21 in the expression cohort,
  302 in the genotype cohort. The 2 missing-genotype subjects are
  expression-only and are silently dropped from the within-T21 regression
  step. Phrase paper text accordingly ("302 of 304").
- **lfcSE column**: chr21 combined output (`deseq2_chr21_genes_both_analyses.csv`)
  does not carry `lfcSE`; the all-genes ploidy-normalized output does.
  The archived forest-plot script joins lfcSE in from the all-genes table.
