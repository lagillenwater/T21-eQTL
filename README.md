# T21-eQTL: Dosage Compensation in Down Syndrome

Cohort-scale extension of Hunter et al. (2023) "Transcription dosage
compensation does not occur in Down syndrome" (BMC Biology 21:228), applied
to the Human Trisome Project (HTP): 304 T21 + 95 Control whole-blood RNA-seq
samples, 302 of the T21 subjects with paired chr21 genotypes.

**Question.** Of the chromosome 21 protein-coding genes that show ploidy-
corrected expression deviations in T21 vs Control, how many can be explained
by common cis-acting allelic variation (eQTLs in GTEx whole blood)?

## Quick start

From the repository root, with R >= 4.2:

```bash
Rscript install_packages.R                        # one-time package install

Rscript scripts/00_preprocess_data.R              # long -> wide count matrix
Rscript scripts/01_deseq2_analysis.R              # trisomy-aware DESeq2
Rscript scripts/02_filter_genotypes.R             # deviating genes + genotype universe
Rscript scripts/03_t21_dosage_boxplots.R          # within-T21 fits + permutation test
Rscript scripts/04_chr21_lane_assignment.R        # MAIN: per-gene lane table
Rscript scripts/05_alluvial_lane_assignment.R     # alluvial flow figure
Rscript scripts/06_chr21_distribution_panel.R     # chr21 vs genome distributions
Rscript scripts/07_three_panel_figure.R           # volcano summary panel
```

Total runtime end-to-end on a laptop: ~30 minutes, dominated by 02 (genotype
streaming) and 03 (per-variant within-T21 regressions). Script 04 also
emits the SankeyMATIC text input - paste into https://sankeymatic.com/build/
to render the lane-flow diagram (instructions in the script's header).

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

## Repository layout

```
T21-eQTL/
  README.md                  # this file (canonical doc)
  install_packages.R
  environment.yml            # conda alternative
  scripts/                   # production pipeline (see Pipeline section)
    lib/                     # shared helpers
    archive/                 # legacy + supplementary scripts (see docs/decisions.md)
  data/                      # inputs - mostly .gitignored
    HTP_WholeBlood_RNAseq_Counts_Synapse.txt   # 3.9 GB raw counts
    P4C_metadata_021921_Costello.txt           # sample metadata
    P4C_Comorbidity_020921.tsv                 # optional comorbidities
    chr21_ds_PASS.csv                          # 6.0 GB T21 chr21 genotypes
    chr21_ctrl_PASS.csv                        # 103 MB Control chr21 genotypes
    GTEx_Analysis_v10_..._Whole_Blood.v10.allpairs.chr21.parquet
                                               # GTEx allpairs (chr21)
    processed/                                 # script outputs - .gitignored
  results/
    tables/                  # CSVs from the pipeline (.gitignored)
    figures/                 # PDFs / PNGs (.gitignored)
  docs/
    decisions.md             # decision log, legacy notes, gotchas
    package_installation_info.txt
```

## Input data

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

Chr21 extract from GTEx v10 allpairs (every cis-window variant tested per
gene, regardless of significance). Filename:
`data/GTEx_Analysis_v10_QTLs_GTEx_Analysis_v10_eQTL_all_associations_Whole_Blood.v10.allpairs.chr21.parquet`.
Script 02 applies `pval_nominal <= 1e-4` to keep the variant universe
manageable.

To re-download the chr21 allpairs, run `bash download_gtex.sh` from the
repo root. It fetches the genome-wide Whole_Blood allpairs file from GTEx
v10 and uses python + pyarrow to write the chr21-only parquet under the
filename above.

## Pipeline

The production pipeline is the 00-07 chain in `scripts/`, backed by shared
helpers in `scripts/lib/` (`cohort.R` - analysis-cohort definition;
`chr21_threshold.R` - chr21-internal robust outlier test, annotation only;
`composition.R` - co-expression composition control; `lane_rules.R` -
the sig_lane classification rule (`assign_sig_lane`); `eqtl_fit.R` -
vectorized per-variant regressions and the gene-level permutation test;
`table1.R` - cohort-characteristics table helpers).

| Script | Purpose | Output |
|---|---|---|
| 00_preprocess_data | Long -> wide gene x sample matrix; match metadata | `data/processed/count_matrix.csv`, `sample_metadata.csv`, `gene_annotations.csv` |
| 01_deseq2_analysis | Trisomy-aware DESeq2 (ploidy normalization matrix; chr21 excluded from size factors; betaPrior=FALSE) | `results/tables/deseq2_all_genes_ploidy_normalized.csv`, `deseq2_chr21_genes_both_analyses.csv`, `deseq2_all_genes_both_analyses.csv`, QC PDFs |
| 02_filter_genotypes | Select deviating chr21 genes by Hunter et al.'s rule (`norm_padj < ALPHA_DE` AND `abs(norm_log2FC) >= DEVIATION_LFC`), restrict to protein-coding, compute the chr21-internal outlier annotation (`dev_z`, `q_outlier`), pull GTEx allpairs cis variants, stream PASS files for those positions | `data/processed/eqtl_supported_genes.csv`, `eqtl_target_variants.csv`, `genotypes_filtered.csv` |
| 03_t21_dosage_boxplots | Per-(variant, gene) within-T21 expression ~ dosage regressions; gene-level cis-eQTL permutation test (`scripts/lib/eqtl_fit.R`); negative-control diagnostics | `results/tables/t21_dosage_per_variant.csv`, `t21_representative_variants.csv`, `eqtl_gene_level_perm.csv`, `eqtl_negative_controls.csv` |
| **04_chr21_lane_assignment** | **Per-gene lane assignment.** Hunter's padj + 1.5-fold rule is the classification split (`sig_lane`, applied by `scripts/lib/lane_rules.R`); deviating genes get a composition control (`scripts/lib/composition.R`) and an `eqtl_lane` terminal from the gene-level permutation test | `results/tables/chr21_lane_assignments.csv`, `chr21_lane_summary.csv`, `chr21_composition_control.csv`, `chr21_k_sensitivity.csv` |
| 05_alluvial_lane_assignment | Alluvial flow (Classification -> Sub-category -> eQTL terminal) + SankeyMATIC export | `results/figures/chr21_lane_alluvial.{pdf,png}`, `results/tables/chr21_lane_sankeymatic_input.txt`, `chr21_lane_alluvial_flow.csv` |
| 06_chr21_distribution_panel | Density + ECDF of chr21 vs baseMean-matched non-chr21 protein-coding distributions; per-lane magnitude scatter | `results/figures/chr21_vs_genome_distribution.{pdf,png}` |
| 07_three_panel_figure | 2x2 volcanoes: A/B all genes before/after ploidy correction, C/D chr21 only; labels read from the lane table | `results/figures/Chr21_DEG.{pdf,png}` |

Legacy and supplementary scripts live under `scripts/archive/`; the
catalog is in [docs/decisions.md](docs/decisions.md).

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



## Outputs

### Tables

`results/tables/`:
- `deseq2_chr21_genes_both_analyses.csv` - chr21 DESeq2 results (raw + ploidy-
  corrected in one row per gene).
- `deseq2_all_genes_ploidy_normalized.csv` - genome-wide ploidy-corrected
  DESeq2 results (reference; classification uses the chr21 table directly).
- **`chr21_lane_assignments.csv`** - canonical per-gene lane table (read
  this for the headline numbers).
- `chr21_lane_summary.csv` - lane counts (all chr21 + after paper filters).
- `chr21_lane_alluvial_flow.csv` - long-format flow data for the alluvial.
- **`chr21_lane_sankeymatic_input.txt`** - SankeyMATIC paste-ready export.
- `t21_dosage_per_variant.csv` - per-variant within-T21 regression fits.

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
- `Chr21_DEG.{pdf,png}` - volcano summary panel (script 07).


## Further documentation

[docs/decisions.md](docs/decisions.md) - decision log (retired
classification filters, the replaced eQTL rule, composition-control null
design), legacy data sources and terminology, the `scripts/archive/`
catalog, and practical gotchas.

## AI assistance

This project utilized the AI assistant Claude, developed by Anthropic, during
the development process. Its assistance included generating initial code
snippets and improving documentation. All AI-generated content was reviewed,
tested, and validated by human developers.

## Citation

Hunter, S., Hendrix, J., Freeman, J., Dowell, R.D., & Allen, M.A. (2023).
Transcription dosage compensation does not occur in Down syndrome.
*BMC Biology* 21:228. https://doi.org/10.1186/s12915-023-01700-4

Original analysis code: https://github.com/Dowell-Lab/DS_Normalization

## License

BSD 3-Clause License — see [LICENSE](LICENSE).
