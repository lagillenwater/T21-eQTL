# T21-eQTL Analysis

## Project Overview

This repository extends the analysis from Hunter et al. (2023) "Transcription
dosage compensation does not occur in Down syndrome" (BMC Biology 21:228) to
the Human Trisome Project (HTP) cohort: 304 T21 + 95 Control whole-blood
RNA-seq samples, 302 of the T21 subjects with paired chr21 genotypes.

**Question.** Of the chromosome 21 protein-coding genes that show ploidy-
corrected expression deviations in T21 vs Control, how many can be explained
by common cis-acting allelic variation (eQTLs in GTEx whole blood) versus
remaining as candidates for true regulatory dosage compensation?

**Headline result** (current pipeline, May 2026):

- 160 chr21 protein-coding genes tested.
- 119 (74%) classified as **Expected dosage** (deviation within 1 cohort-noise SD).
- 41 outside cohort noise:
  - 6 DE_high, all eQTL-supported (TSPEAR, MX1, RIPK4, COL6A2, CYYR1, YBEY).
  - 15 DE_low: 9 eQTL-supported, 3 eQTL-tested-not-supported (BACE2, ADARB1,
    NRIP1), 3 with no GTEx eQTL coverage (LTN1, RUNX1, ZBTB21).
  - 20 in High repeats / Low expression / Not DE outside cohort noise.

---

## Quick Start

From the repository root, with R >= 4.2:

```bash
Rscript install_packages.R                        # one-time package install

Rscript scripts/00_preprocess_data.R              # long -> wide count matrix
Rscript scripts/01_deseq2_analysis.R              # trisomy-aware DESeq2
Rscript scripts/02_filter_genotypes.R             # deviating-gene selection (Hunter's padj + 1.5-fold rule) + genotype universe
Rscript scripts/03_t21_dosage_boxplots.R          # within-T21 per-variant fits + gene-level eQTL permutation test
Rscript scripts/04_chr21_lane_assignment.R        # MAIN: per-gene lane table + composition control + SankeyMATIC export
Rscript scripts/05_alluvial_lane_assignment.R     # alluvial flow figure + SankeyMATIC export
Rscript scripts/06_chr21_distribution_panel.R     # chr21 vs genome distributions
Rscript scripts/07_three_panel_figure.R           # 3-panel summary figure
```



Total runtime end-to-end on a laptop: ~30 minutes, dominated by 02 (genotype
streaming) and 03 (per-variant within-T21 regressions). Script 04 also
emits the SankeyMATIC text input - paste into https://sankeymatic.com/build/
to render the diagram (instructions in the script's header).

---

## Installation

### Method 1: Automated (recommended)

```bash
Rscript install_packages.R
```

Installs CRAN + Bioconductor packages, writes
`docs/package_installation_info.txt`. ~5-15 minutes.

### Method 2: Manual

```r
install.packages(c(
  "tidyverse", "data.table", "ggplot2", "ggrepel",
  "ggalluvial", "patchwork", "RColorBrewer", "viridis",
  "here", "arrow"
))
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("DESeq2")
```

### Compilation toolchain (if a source install fails)

- macOS: `xcode-select --install`
- Linux (Debian/Ubuntu): `sudo apt-get install r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev`
- Windows: install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)

### Verification

```r
for (p in c("tidyverse", "data.table", "DESeq2", "ggalluvial",
            "patchwork", "arrow")) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("%-15s %s\n", p, if (ok) "OK" else "FAILED"))
}
```

Tested versions: R 4.5.2; tidyverse 2.0; data.table 1.16; DESeq2 1.46;
ggalluvial 0.12; arrow 17.

---

## Repository Layout

```
T21-eQTL/
  CLAUDE.md                  # this file (canonical doc)
  README.md                  # slim entry point
  install_packages.R
  environment.yml            # conda alternative
  scripts/                   # production pipeline (see Pipeline section)
  data/                      # inputs - mostly .gitignored
    HTP_WholeBlood_RNAseq_Counts_Synapse.txt   # 3.9 GB raw counts
    P4C_metadata_021921_Costello.txt           # sample metadata
    P4C_Comorbidity_020921.tsv                 # optional comorbidities
    chr21_ds_PASS.csv                          # 6.0 GB T21 chr21 genotypes
    chr21_ctrl_PASS.csv                        # 103 MB Control chr21 genotypes
    Whole_Blood.v10.eQTLs.signif_pairs.parquet # legacy GTEx significant pairs
    GTEx_Analysis_v10_..._Whole_Blood.v10.allpairs.chr21.parquet
                                               # current GTEx allpairs (chr21)
    processed/                                 # script outputs - .gitignored
  results/
    tables/                  # CSVs from the pipeline (.gitignored)
    figures/                 # PDFs / PNGs (.gitignored)
  docs/
    package_installation_info.txt
```

---

## Input Data

### 1. RNA-seq counts (HTP whole blood)

**File**: `data/HTP_WholeBlood_RNAseq_Counts_Synapse.txt` (3.9 GB, long format)

Columns: `LabID`, `Sample_type`, `Platform`, `EnsemblID`, `Gene_name`, `Chr`,
`Gene_type`, `Units`, `Value`. ~24M rows = ~400 samples x ~60k genes.
Script 00 pivots this to a gene x sample matrix (`data/processed/count_matrix.csv`).

### 2. Sample metadata

**File**: `data/P4C_metadata_021921_Costello.txt` (587 samples). Columns include
`RecordID`, `Sex`, `Karyotype` (T21 / Control / other), `LabID` (without the
tissue suffix used in the count file), `Age_at_visit`, `BMI`,
`Sample_source`. Matched to count data by stripping the trailing letter
suffix. After matching: 304 T21 + 95 Control with both expression and
metadata.

### 3. Comorbidity data (optional)

**File**: `data/P4C_Comorbidity_020921.tsv`. Not currently used downstream;
kept for future covariate adjustment.

### 4. HTP chr21 genotypes

**Files**: `data/chr21_ds_PASS.csv` and `data/chr21_ctrl_PASS.csv`. VCF-style
CSVs filtered to PASS variants, chr21 only. T21 file uses ploidy-3 calls;
Control file uses ploidy-2. Streamed by script 02 to extract only the
positions of the eQTL universe.

### 5. GTEx whole-blood eQTLs

**Current source**: chr21 extract from GTEx v10 allpairs (every cis-window
variant tested per gene, regardless of significance). Filename:
`data/GTEx_Analysis_v10_QTLs_GTEx_Analysis_v10_eQTL_all_associations_Whole_Blood.v10.allpairs.chr21.parquet`.

**Legacy source**: `data/Whole_Blood.v10.eQTLs.signif_pairs.parquet`
(per-gene FDR-passing pairs only). Older runs of script 02 used this; the
current pipeline reads allpairs and applies `pval_nominal <= 1e-4` to keep
the variant universe manageable.

To re-download the chr21 allpairs, run `bash download_gtex.sh` from the
repo root. It fetches the genome-wide Whole_Blood allpairs file from GTEx
v10 and uses python + pyarrow to write the chr21-only parquet under the
filename above.

---

## Pipeline

The production pipeline is the 00-07 chain in `scripts/`, backed by shared
helpers in `scripts/lib/` (`cohort.R` - analysis-cohort definition;
`chr21_threshold.R` - chr21-internal robust outlier test, annotation only;
`composition.R` - co-expression composition control; `lane_rules.R` -
the sig_lane classification rule (`assign_sig_lane`); `eqtl_fit.R` -
vectorized per-variant regressions and the gene-level permutation test;
`table1.R` - cohort-characteristics table helpers).

### Production pipeline (scripts/00-07)

| Script | Purpose | Output |
|---|---|---|
| 00_preprocess_data | Long -> wide gene x sample matrix; match metadata | `data/processed/count_matrix.csv`, `sample_metadata.csv`, `gene_annotations.csv` |
| 01_deseq2_analysis | Trisomy-aware DESeq2 (ploidy normalization matrix; chr21 excluded from size factors; betaPrior=FALSE) | `results/tables/deseq2_all_genes_ploidy_normalized.csv`, `deseq2_chr21_genes_both_analyses.csv`, `deseq2_all_genes_both_analyses.csv`, QC PDFs |
| 02_filter_genotypes | Select deviating chr21 genes by Hunter et al.'s rule (`norm_padj < ALPHA_DE` AND `abs(norm_log2FC) >= DEVIATION_LFC`), restrict to protein-coding, compute the chr21-internal FDR-outlier annotation (`dev_z`, `q_outlier`, via `scripts/lib/chr21_threshold.R`), pull GTEx allpairs cis variants (`pval_nominal <= 1e-4`), stream PASS files for those positions | `data/processed/eqtl_supported_genes.csv`, `eqtl_target_variants.csv`, `genotypes_filtered.csv` |
| 03_t21_dosage_boxplots | Per-(variant, gene) within-T21 expression ~ dosage regressions; gene-level cis-eQTL permutation test (`scripts/lib/eqtl_fit.R`); negative controls (direction-flip, genotype-permutation) on the retired any-variant rule | `results/tables/t21_dosage_per_variant.csv`, `t21_representative_variants.csv`, `eqtl_gene_level_perm.csv`, `eqtl_negative_controls.csv` |
| **04_chr21_lane_assignment** | **Per-gene lane assignment.** Hunter's padj + 1.5-fold rule is the classification split (`sig_lane`: `Expected_dosage`/`High_repeats`/`Low_expression`/`DE_low`/`DE_high`, applied by `scripts/lib/lane_rules.R`); deviating genes get a composition control (`scripts/lib/composition.R`, PROGRAM/MIXED/GENE-SPECIFIC verdict) and an `eqtl_lane` terminal from the gene-level permutation test (`cis_eqtl`/`no_cis_eqtl`/`no_GTEx_data`) | `results/tables/chr21_lane_assignments.csv`, `chr21_lane_summary.csv`, `chr21_composition_control.csv`, `chr21_k_sensitivity.csv` |
| 05_alluvial_lane_assignment | Alluvial flow (Classification -> Sub-category -> eQTL terminal) + SankeyMATIC export | `results/figures/chr21_lane_alluvial.{pdf,png}`, `results/tables/chr21_lane_sankeymatic_input.txt`, `chr21_lane_alluvial_flow.csv` |
| 06_chr21_distribution_panel | Density + ECDF of chr21 vs baseMean-matched non-chr21 protein-coding distributions; per-lane magnitude scatter | `results/figures/chr21_vs_genome_distribution.{pdf,png}` |
| 07_three_panel_figure | 2x2 volcanoes: A/B all genes before/after ploidy correction, C/D chr21 only; labels read from the lane table. Lane flow rendered separately via SankeyMATIC from chr21_lane_sankeymatic_input.txt | `results/figures/Chr21_DEG.{pdf,png}` |

### Archived / supplementary (scripts/archive/, do not edit)

Original paper-style Panel D pipeline:
- `02_categorize_genes` - paper-style gene categorization
- `03_volcano_plot` - diagnostic volcano
- `04_alluvial_plot` - original Panel D Sankey
- `05_eqtl_analysis` - original eQTL cross-reference (template-based)
- `06_alluvial_with_eqtl` - enhanced Panel D with eQTL terminals
- `07_sankeymatic_export` - SankeyMATIC export from old categorization
- `08_expected_dosage_eqtl` - eQTL analysis for >=1.5 FC genes

Supplementary outputs from the cohort-scale chain:
- `10_eqtl_genotype_concordance` - per-variant T21 vs Control dosage
   means; directional concordance (`results/tables/eqtl_genotype_concordance_*.csv`)
- `14_dosage_lane_boxplots` - per-quadrant focused boxplot PDFs
- `16_chr21_de_forest_plot` - per-DE-gene log2FC + 95% CI
- `17_chr21_quadrant_plot` - within-T21 slope vs deviation scatter

Exploratory / one-offs:
- `diagnostic_check`, `investigate_pc2`, `pca_chr21_only`,
  `process_blacklist`, `run_all.sh` (the old shell driver)

---

## Methodology

### Trisomy-aware DESeq2

The single biggest methodological correction from the paper. Standard DESeq2
(default null `FC = 1`) calls every chr21 gene differentially expressed in
T21 because the expected FC is `1.5`, not `1`. Two issues compound:

1. **Wrong null hypothesis.** Significance tests against `FC = 1` over-call DE
   on chr21 even when expression is exactly proportional to copy number.
2. **Inflated dispersion.** Including chr21 in the size-factor estimation
   shrinks all fold changes toward 1, masking real signal elsewhere.

Fixes (script 01):
- Build a per-gene-per-sample ploidy-normalization matrix: 1.5 for
  chr21 genes in T21 samples, 1.0 elsewhere. Pass via DESeq2's `normMatrix`.
- Exclude chr21 from size-factor computation.
- `betaPrior = FALSE` (no shrinkage), so the MAP fold change reflects the
  raw maximum-likelihood estimate.
- Apply a minimum baseMean filter (paper's "second quintile"; pipeline uses
  the 20th-percentile baseMean cutoff for consistency with script 02).

After ploidy normalization, the appropriate null on chr21 is back to `FC = 1`,
so DESeq2's standard p-value testing applies cleanly.

### Hunter's classification rule

The current classification rule matches Hunter et al.'s own criterion
directly, rather than a cohort-derived noise threshold: a chr21 gene
deviates if `norm_padj < ALPHA_DE` (0.01) AND `abs(norm_log2FC) >=
DEVIATION_LFC` (`log2(1.5)`, their FC >= 1.5 cut, applied on the
ploidy-corrected scale). That is tier 1, the primary result. A labelled
secondary tier 2 (`DEVIATION_LFC_T2 = log2(4/3)`) admits near-threshold
genes with the same padj requirement; both tiers get the composition
control and the cis-eQTL permutation test, and the `tier` column in
`chr21_lane_assignments.csv` records which is which. This is the sole gate on `sig_lane` /
`passes_magnitude_filter` in scripts 02 and 04.

A chr21-internal FDR-controlled robust outlier test (`scripts/lib/
chr21_threshold.R`: median/MAD null estimated from expressed, non-repeat
chr21 genes, BH-FDR on the robust z-score) is computed alongside and
written as annotation columns `dev_z` / `q_outlier` (plus
`chr21_k_sensitivity.csv`), but does not gate classification. Ploidy
normalization only acts on chr21 genes, so a non-chr21 cohort-noise SD is
not a valid reference for this null; retiring the cohort-SD filter removed
that mismatch.

### Composition control

`scripts/lib/composition.R`, run in script 04 for every deviating gene.
A chr21 gene can appear to deviate because the blood cell population that
expresses it shifted in T21, not because of a regulatory effect on the
gene itself - a composition shift moves a whole co-expression program, not
one gene. For each deviating gene: find its 20 most-correlated non-chr21
partners (correlated in controls only, so the karyotype effect cannot leak
into the neighborhood definition), take their median T21-vs-Control
log2FC (`partner_lfc`), and compare against a **correlation-matched** null
to get `p_partners`.

The null must be matched, not independent. `partner_null(L_ctrl, lfc,
n_partners = 20, n_draw = 300, seed = 1)` draws `n_draw` random seed genes
from the same non-chr21 expressed pool and, for each seed, takes the median
log2FC of *its own* top-`n_partners` correlated genes - a null draw is a
co-expression module built exactly like the observed one. The seed is
excluded from its own partner set. An earlier version drew independent
random 20-gene sets, which is anti-conservative: a module moves together,
so its median log2FC is several times more variable (measured here: SD
0.259 matched vs 0.058 independent), and an independent null called almost
any partner shift significant. `partner_p()` reports the one-sided
empirical p as `(1 + k) / (n_draw + 1)`, so it is never exactly 0.
Verdict: **PROGRAM** if
`p_partners < 0.05` and `partner_lfc / gene_lfc >= 0.5`; **MIXED** if
`p_partners < 0.05` only; else **GENE-SPECIFIC**. Written to
`results/tables/chr21_composition_control.csv`
(`Gene_name, gene_lfc, partner_lfc, p_partners, program_share,
residual_lfc, verdict`), with `residual_lfc = gene_lfc - partner_lfc`
merged back into `chr21_lane_assignments.csv`.

### Lane assignment

`scripts/04_chr21_lane_assignment.R` produces `chr21_lane_assignments.csv` -
the canonical per-gene table. Each row has:

- DESeq2 stats: `baseMean`, `raw_log2FC`, `raw_FC`, `norm_log2FC`, `norm_padj`.
- Magnitude annotation: `deviation_magnitude`, `eligible_idx` (not
  repeat-flagged, not low-expression, corrected log2FC estimable),
  `passes_magnitude_filter` (`eligible_idx` AND Hunter's cut), plus the
  annotation-only outlier columns `dev_z`, `q_outlier`.
- Categorization flags: `low_expr`, `high_repeat`.
- **`sig_lane`**: one of `Expected_dosage` (fails Hunter's rule),
  `High_repeats`, `Low_expression`, `DE_low`, `DE_high`, `Not_DE_outside_noise`.
- Composition control: `verdict` (PROGRAM/MIXED/GENE-SPECIFIC), `residual_lfc`.
- **`eqtl_lane`**: one of `not_evaluated` (non-deviating), `cis_eqtl`
  (gene-level permutation `q_gene_bh < FDR_GENE`), `no_cis_eqtl`,
  `no_GTEx_data`. A detection result, not an "explained by eQTL" claim -
  see `docs/REPO_STATE.md` decision log.
- Gene-level permutation result: `p_gene_perm`, `q_gene_bh`,
  `cis_eqtl_detected`, `best_variant`.
- Locus-level eQTL aggregation (retained for context, not the
  classification rule): `n_cis_total`, `n_dir_match`, `n_supp_with_repro`.
- Representative-variant choices for downstream visualization:
  `strongest_supp_variant`, `strongest_dir_variant`,
  `strongest_overall_variant`.

### eQTL detection: gene-level permutation test

The locus-level "any cis variant matches direction and reproduces in T21"
rule (script 03's `strongest_supp_variant` logic) is retained as context
but is **not** the classification rule: `results/tables/
eqtl_negative_controls.csv` shows it returns 100% "explained" even when
the observed deviation direction is artificially flipped, so it does not
discriminate real signal from chance. The classification rule is a
gene-level permutation test (`scripts/lib/eqtl_fit.R`, run in script 03):
for each deviating gene, the best-variant test statistic is compared
against its null distribution under permutation of genotype-to-expression
assignment, giving `p_gene_perm`; BH-adjusted across deviating genes to
`q_gene_bh`. `eqtl_lane = cis_eqtl` when `q_gene_bh < FDR_GENE` (0.05).

### Constants worth knowing

Defined at the top of scripts 02, 03, and 04; change once, propagates through:

- `ALPHA_DE = 0.01` - paper's padj threshold (Fig. 2B, 3B).
- `DEVIATION_LFC = log2(1.5)` - Hunter's FC >= 1.5 cut on the
  ploidy-corrected scale; the sole classification threshold. The one dial
  in this analysis - see `docs/REPO_STATE.md` decision log.
- `ALPHA_REPRO = 0.05` - within-T21 nominal p for a cis variant to count as
  reproducible (locus-level context column, not the classification rule).
- `FDR_GENE = 0.05` - BH threshold on `q_gene_bh` for `cis_eqtl_detected`.
- `OUTLIER_FDR = 0.10` - FDR threshold for the chr21-internal outlier
  annotation (`dev_z`/`q_outlier`); **annotation only, retired from
  classification** (was `MAGNITUDE_THRESHOLD` in an earlier cohort-SD
  filter, also retired).
- `LOW_EXPR_QUANT = 0.20` - baseMean q20 (paper's "second quintile") low-
  expression filter.
- `GTEX_PVAL_KEEP = 1e-4` - nominal cis-eQTL pval cutoff in allpairs
  (matches the effective signif_pairs cutoff).
- `RESTRICT_TO_PROTEIN_CODING = TRUE` - restrict the target chr21 gene set
  (and the eligible-gene pool for the chr21-internal outlier null) to
  protein-coding genes.

---

## Outputs

### Tables

`results/tables/`:
- `deseq2_chr21_genes_both_analyses.csv` - chr21 DESeq2 results (raw + ploidy-
  corrected in one row per gene).
- `deseq2_all_genes_ploidy_normalized.csv` - genome-wide ploidy-corrected
  DESeq2 results (all genes, used for the ploidy-corrected volcano panel and
  general reference; not used for gene classification, which is Hunter's
  padj + 1.5-fold rule on the chr21 table directly).
- **`chr21_lane_assignments.csv`** - canonical per-gene lane table (read
  this for the headline numbers).
- `chr21_lane_summary.csv` - lane counts (all chr21 + after paper filters).
- `chr21_lane_alluvial_flow.csv` - long-format flow data for the alluvial.
- **`chr21_lane_sankeymatic_input.txt`** - SankeyMATIC paste-ready export.
- `t21_dosage_per_variant.csv` - per-variant within-T21 regression fits.
- `t21_representative_variants.csv` - per-gene strongest supportive variant
  (legacy, used by `scripts/archive/14_dosage_lane_boxplots.R`).

`data/processed/`:
- `count_matrix.csv` - gene x sample expression matrix.
- `sample_metadata.csv` - matched, filtered metadata.
- `gene_annotations.csv` - chr / gene_type lookup.
- `eqtl_supported_genes.csv` - target gene roster for eQTL stage.
- `eqtl_target_variants.csv` - cis variants per target gene (from GTEx).
- `genotypes_filtered.csv` - HTP genotypes at the cis variants.

### Figures

`results/figures/`:
- **`chr21_lane_alluvial.{pdf,png}`** - main lane-flow visualization
  (script 05).
- `chr21_vs_genome_distribution.{pdf,png}` - chr21 vs non-chr21
  ploidy-corrected log2FC distributions + per-lane magnitude scatter
  (script 06).
- Supplementary figures from `scripts/archive/` (run them on demand):
  per-quadrant boxplots (`chr21_dosage_de_*.pdf`), DE forest plot
  (`chr21_de_forest_plot.{pdf,png}`), within-T21 quadrant plot
  (`chr21_quadrant_plot.{pdf,png}`), plus the paper-style legacy outputs.

---

## Running and re-running

The pipeline is idempotent: running it from scratch reproduces the same
outputs. Stages are independent at the file level - if you change script
04, just re-run 04 -> 05 -> 06 (not 00 -> 03).

Common re-run patterns:

- **Tweak the deviation threshold**: edit `DEVIATION_LFC` at top of
  scripts 02 and 04, then re-run 02 -> 03 -> 04 -> 05 -> 06 -> 07. (Cheaper
  if only the lane assignment is changing: re-run 04 -> 05 -> 06 -> 07 only.)
- **Switch to/from protein-coding**: set `RESTRICT_TO_PROTEIN_CODING` in
  scripts 02, 04, 06. Re-run from 02.
- **Change padj threshold**: `ALPHA_DE` in scripts 02 and 04. Re-run 02 onward.

The cost concentrations: script 02 (~3 min, awk-streaming PASS files);
script 03 (~5-15 min depending on variant count, per-variant regressions).
The visualization scripts (05, 06) finish in seconds each.

---

## Common gotchas

- **`Sig_high_FC` vs `DE_high`**: an early version of the pipeline used
  `Sig_high_FC` in `eqtl_supported_genes.csv` `gene_set` column (written
  by script 02). The current lane terminology (in script 04 onward) is
  `DE_high`. Both refer to the same set: genes with `norm_padj < 0.01` AND
  `norm_log2FC >= log2(1.5)`. The cut is on the PLOIDY-CORRECTED log2FC and
  split by its sign, not on `raw_FC` - the raw chr21 fold change centres on
  1.5 by ploidy alone, so a raw-FC cut would select on the trisomy itself.
  Figure labels that said "raw FC" were wrong and were corrected.
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

### Failed-experiment scripts

`scripts/archive/diagnostic_check.R`, `investigate_pc2.R`,
`pca_chr21_only.R`, and `process_blacklist.R` are exploratory / one-off
scripts. Kept but not expected to run cleanly against the current state.

---

## Paper reference

Hunter, S., Hendrix, J., Freeman, J., Dowell, R.D., & Allen, M.A. (2023).
Transcription dosage compensation does not occur in Down syndrome.
*BMC Biology* 21:228. https://doi.org/10.1186/s12915-023-01700-4

Original analysis code: https://github.com/Dowell-Lab/DS_Normalization
