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
Rscript scripts/02_filter_genotypes.R             # gene + variant + genotype universe
Rscript scripts/03_t21_dosage_boxplots.R          # within-T21 per-variant fits
Rscript scripts/04_chr21_lane_assignment.R        # MAIN: per-gene lane table + SankeyMATIC export
Rscript scripts/05_chr21_distribution_panel.R     # chr21 vs genome distributions
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
Control file uses ploidy-2. Streamed by script 09 to extract only the
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

The production pipeline is the 00-05 chain in `scripts/`.

### Production pipeline (scripts/00-05)

| Script | Purpose | Output |
|---|---|---|
| 00_preprocess_data | Long -> wide gene x sample matrix; match metadata | `data/processed/count_matrix.csv`, `sample_metadata.csv`, `gene_annotations.csv` |
| 01_deseq2_analysis | Trisomy-aware DESeq2 (ploidy normalization matrix; chr21 excluded from size factors; betaPrior=FALSE) | `results/tables/deseq2_chr21_combined.csv`, `deseq2_all_genes_ploidy_normalized.csv`, QC PDFs |
| 02_filter_genotypes | Apply magnitude filter to gene set (1.0 cohort-SD threshold), restrict to protein-coding, pull GTEx allpairs cis variants (`pval_nominal <= 1e-4`), stream PASS files for those positions | `data/processed/eqtl_supported_genes.csv`, `eqtl_target_variants.csv`, `genotypes_filtered.csv` |
| 03_t21_dosage_boxplots | Per-(variant, gene) within-T21 expression ~ dosage regressions | `results/tables/t21_dosage_per_variant.csv`, `t21_representative_variants.csv` |
| **04_chr21_lane_assignment** | **Per-gene lane assignment.** Cohort-noise filter is the first split; survivors get DE_low/DE_high/High_repeats/Low_expression/Not_DE classification; DE genes get eQTL-supported / eQTL-tested-not-supported / no_GTEx_data terminal | `results/tables/chr21_lane_assignments.csv`, `chr21_lane_summary.csv` |
| 05_alluvial_lane_assignment | Alluvial flow (Cohort filter -> Sub-category -> eQTL terminal) + SankeyMATIC export | `results/figures/chr21_lane_alluvial.{pdf,png}`, `results/tables/chr21_lane_sankeymatic_input.txt`, `chr21_lane_alluvial_flow.csv` |
| 06_chr21_distribution_panel | Density + ECDF of chr21 vs baseMean-matched non-chr21 protein-coding distributions; per-lane magnitude scatter | `results/figures/chr21_vs_genome_distribution.{pdf,png}` |

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
  the 20th-percentile baseMean cutoff for consistency with script 09).

After ploidy normalization, the appropriate null on chr21 is back to `FC = 1`,
so DESeq2's standard p-value testing applies cleanly.

### Cohort-noise (within-cohort SD) filter

The paper's family-of-4 dataset got an effective magnitude filter for free
(small `n` -> big lfcSE -> only large effects survive `padj < 0.01`). At our
n = 304 + 95, padj catches arbitrarily small deviations as significant -
including ones smaller than typical cohort variation.

The filter (constant `MAGNITUDE_THRESHOLD = 1.0` in scripts 02 and 04):

1. From `deseq2_all_genes_ploidy_normalized.csv`, take all non-chr21
   protein-coding genes. After ploidy correction these are expected to
   center near `log2FC = 0`; their SD is the cohort-noise yardstick
   (currently `0.40` in log2 units).
2. For each chr21 gene, compute `|norm_log2FC| / cohort_sd`.
3. Genes with ratio < 1.0 are categorized as **Expected dosage** before any
   eQTL evaluation - their statistically significant deviations are within
   typical genome-wide cohort variation, so no eQTL story is sought.
4. Genes with ratio >= 1.0 (currently 41 of 160 protein-coding chr21 genes)
   proceed to further categorization.

### Lane assignment

`scripts/04_chr21_lane_assignment.R` produces `chr21_lane_assignments.csv` -
the canonical per-gene table. Each row has:

- DESeq2 stats: `baseMean`, `raw_log2FC`, `raw_FC`, `norm_log2FC`, `norm_padj`.
- Magnitude annotation: `deviation_magnitude`, `deviation_vs_cohort_sd`,
  `passes_magnitude_filter`.
- Categorization flags: `low_expr`, `high_repeat`.
- **`sig_lane`**: one of `Expected_dosage` (cohort-noise filter failure),
  `High_repeats`, `Low_expression`, `DE_low`, `DE_high`, `Not_DE_outside_noise`.
- **`eqtl_lane`**: one of `not_evaluated` (non-DE), `explained` (>=1 cis
  variant matches deviation direction in GTEx AND reproduces in T21 at
  `t21_p < 0.05`), `unexplained`, `no_GTEx_data`.
- Locus-level eQTL aggregation: `n_cis_total`, `n_dir_match`,
  `n_supp_with_repro`.
- Three representative-variant choices for downstream visualization:
  `strongest_supp_variant` (T21-reproducible direction match),
  `strongest_dir_variant` (best GTEx direction-match candidate, used when
  `strongest_supp_variant` is missing), `strongest_overall_variant`
  (smallest GTEx pval regardless of direction).

### eQTL "supported" definition

A cis variant counts as supportive of a gene's deviation if both:
1. `sign(GTEx_slope) == sign(norm_log2FC)` (the eQTL points the same way as
   the observed deviation - paper's directional explanation criterion).
2. `sign(within_T21_slope) == sign(GTEx_slope)` AND `within_T21_p < 0.05`
   (the eQTL reproduces in our T21 cohort at nominal significance, not just
   directional concordance).

A gene is `eqtl_lane = explained` if at least one cis variant in its GTEx
allpairs window satisfies both criteria.

### Constants worth knowing

Defined at the top of scripts 02 and 04; change once, propagates through:

- `ALPHA = 0.01` - paper's padj threshold (Fig. 2B, 3B).
- `ALPHA_REPRO = 0.05` - within-T21 nominal p for a cis variant to count as
  reproducible.
- `MAGNITUDE_THRESHOLD = 1.0` - cohort-noise SD threshold for the magnitude
  filter.
- `LOW_EXPR_QUANT = 0.20` - baseMean q20 (paper's "second quintile") low-
  expression filter.
- `GTEX_PVAL_KEEP = 1e-4` - nominal cis-eQTL pval cutoff in allpairs
  (matches the effective signif_pairs cutoff).
- `RESTRICT_TO_PROTEIN_CODING = TRUE` - restrict both target chr21 set and
  cohort-noise reference to protein-coding genes.

---

## Outputs

### Tables

`results/tables/`:
- `deseq2_chr21_combined.csv` - chr21 DESeq2 results (raw + ploidy-corrected
  in one row per gene).
- `deseq2_all_genes_ploidy_normalized.csv` - genome-wide ploidy-corrected
  DESeq2 results; also the source of the cohort-noise SD reference.
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

- **Tweak the magnitude threshold**: edit `MAGNITUDE_THRESHOLD` at top of
  scripts 02 and 04, then re-run 02 -> 03 -> 04 -> 05 -> 06. (Cheaper if
  only the lane assignment is changing: re-run 04 -> 05 -> 06 only.)
- **Switch to/from protein-coding**: set `RESTRICT_TO_PROTEIN_CODING` in
  scripts 02, 04, 06. Re-run from 02.
- **Change padj threshold**: `ALPHA` in script 04. Re-run 04 onward.

The cost concentrations: script 02 (~3 min, awk-streaming PASS files);
script 03 (~5-15 min depending on variant count, per-variant regressions).
The visualization scripts (05, 06) finish in seconds each.

---

## Common gotchas

- **`Sig_high_FC` vs `DE_high`**: an early version of the pipeline used
  `Sig_high_FC` in `eqtl_supported_genes.csv` `gene_set` column (written
  by script 02). The current lane terminology (in script 04 onward) is
  `DE_high`. Both refer to the same set: genes with `raw_FC >= 1.5` AND
  `norm_padj < 0.01` (after the cohort-SD filter is applied first).
- **Self-loops in SankeyMATIC**: `chr21_lane_sankeymatic_input.txt` skips
  level-to-level passes where the source and target name are identical
  (e.g., the 119 Expected dosage genes terminate at level 2 and would
  otherwise self-loop at levels 3 and 4). The flow file
  `chr21_lane_alluvial_flow.csv` keeps them.
- **Two T21 expression-only subjects**: 304 T21 in the expression cohort,
  302 in the genotype cohort. The 2 missing-genotype subjects are
  expression-only and are silently dropped from the within-T21 regression
  step. Phrase paper text accordingly ("302 of 304").
- **lfcSE column**: chr21 combined output (`deseq2_chr21_combined.csv`)
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
