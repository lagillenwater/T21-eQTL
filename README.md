# T21-eQTL: Dosage Compensation in Down Syndrome

Cohort-scale extension of Hunter et al. (2023) "Transcription dosage
compensation does not occur in Down syndrome" (BMC Biology 21:228), applied
to the Human Trisome Project (HTP): 304 T21 + 95 Control whole-blood RNA-seq
samples, 302 of the T21 subjects with paired chr21 genotypes.

**Question.** Of the chromosome 21 genes (protein-coding, lncRNA, pseudogene) that show ploidy-
corrected expression deviations in T21 vs Control, how many carry a
detectable common cis-eQTL in GTEx whole blood (gene-level permutation
test), and how many do not and so remain open as candidates for regulatory
dosage compensation? A cis-eQTL call is a detection result; it does not by
itself establish that the eQTL accounts for the deviation.

## Quick start

From the repository root, with R >= 4.2:

```bash
Rscript install_packages.R                        # one-time package install

Rscript scripts/00_preprocess_data.R                            # long -> wide count matrix; karyotype subtype
Rscript scripts/01_deseq2_analysis.R --run baseline             # trisomy-aware DESeq2; cohort; expression artifact
Rscript scripts/02_filter_genotypes.R --run baseline            # deviating-gene selection + genotype universe
Rscript scripts/03_t21_dosage_boxplots.R --run baseline         # within-T21 fits + gene-level eQTL permutation test
Rscript scripts/04_chr21_lane_assignment.R --run baseline       # MAIN: per-gene lane table
Rscript scripts/05_sankeymatic_export.R --run baseline          # SankeyMATIC input for the lane-flow figure
Rscript scripts/06_chr21_distribution_panel.R --run baseline    # ploidy-correction effect: chr21 vs chr22
Rscript scripts/07_three_panel_figure.R --run baseline          # volcano figures: all genes, chr21 only
Rscript scripts/11_eqtl_figures.R --run baseline                # eQTL dosage panels, chr21 map
Rscript scripts/13_af_deviation_bound.R --run baseline         # can the cis-eQTLs move the group ratio? (AF bound)
Rscript scripts/15_spike_in_power.R --run baseline             # power per deviating gene by spike-in at its GTEx aFC
Rscript scripts/16_base_rate_expected_dosage.R --run baseline  # the same eQTL test on every Expected-dosage gene (base rate)
Rscript scripts/17_distal_inflation.R --run baseline           # negative control: lambda of distal genotype-expression tests
Rscript scripts/18_detection_vs_effect_size.R --run baseline   # deviating vs Expected-dosage genes at matched effect size; after 15, 16
Rscript scripts/19_gtex_tissue_expression.R --run baseline     # GTEx v10 tissue expression of every chr21 gene, against the eQTLs; after 16
Rscript scripts/S1_covariate_evidence.R --run baseline          # supplement: evidence for covariates and exclusions

# Covariate- and composition-adjusted run: same scripts with --run adjusted, then
Rscript scripts/10_compare_runs.R --run adjusted --against baseline   # per-gene comparison of the two runs
Rscript scripts/audit_baseline_drift.R --run baseline --archive results/archive/2026-09-04_flat
```

Total runtime end-to-end on a laptop: ~30 minutes, dominated by 02 (genotype
streaming) and 03 (per-variant within-T21 regressions). Script 05 also
emits the SankeyMATIC text input
(`results/runs/<name>/tables/chr21_lane_sankeymatic_input.txt`); paste it into
https://sankeymatic.com/build/
to render the lane-flow diagram.

## Runs

Every pipeline script from 01 onward takes `--run <name>` (default
`baseline`; the `T21_RUN` environment variable is honoured when the flag is
absent) and reads `config/runs/<name>.R`: covariates, composition source,
cohort exclusions, ploidy, and every classification threshold. Outputs go
to `results/runs/<name>/{processed,tables,figures}`. Two runs are defined:

- `baseline`: design `~ karyotype`, no exclusions. The primary
  result. Regenerated through the run mechanism and checked against the
  2026-09-04 flat outputs (`results/archive/2026-09-04_flat/`, historical)
  by `scripts/audit_baseline_drift.R`.
- `adjusted`: design `~ age + sex + BMI + sample source + 19 CyTOF cell
  fractions + karyotype`, mosaic T21 subjects excluded. Requires the CyTOF
  table (`download_synapse_celltypes.sh`).

Script 01 writes the run's cohort (`processed/analysis_cohort.csv`, plus
`cohort_roster.csv` with an `excluded_reason` per dropped sample) and one
expression artifact, `processed/expression_adjusted.csv`: log2-CPM
(library size excluding chr21) with the run's covariate effects removed and
karyotype kept (`scripts/lib/covariates.R::adjust_expression`). Script 03
reads that artifact, so the eQTL fits are on adjusted data in the
adjusted run. In baseline nothing is removed. The artifact is not
ploidy-scaled; between-group fold changes come from DESeq2.

To add a covariate, list it in `covariates` (it must be a column of
`data/processed/sample_metadata.csv` named in `RUN_COVARIATES` in
`scripts/lib/run.R`). To exclude samples, add `column = c(values)` to
`exclude`. Thresholds live only in the config. Re-run from script 01.

`scripts/S1_covariate_evidence.R --run <name>` writes the evidence
supplement (covariates vs karyotype; cell-fraction reach into expression;
mosaic subtype vs expression; per-covariate attribution of the fold-change
shift). `scripts/10_compare_runs.R --run adjusted --against baseline`
writes the per-gene run comparison and lane-transition table.

## Installation

### Method 1: Automated (recommended)

```bash
Rscript install_packages.R
```

Installs CRAN + Bioconductor packages, writes
`docs/package_installation_info.txt`. ~5-15 minutes.

### Method 2: Manual

```r
pkgs <- c("tidyverse", "data.table", "ggplot2", "ggrepel",
          "patchwork", "RColorBrewer", "viridis",
          "here", "arrow")
install.packages(pkgs)
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("DESeq2")
```

### Verification

Checks every package from the manual block plus DESeq2:

```r
pkgs <- c("tidyverse", "data.table", "ggplot2", "ggrepel",
          "patchwork", "RColorBrewer", "viridis",
          "here", "arrow", "DESeq2")
for (p in pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("%-15s %s\n", p, if (ok) "OK" else "FAILED"))
}
```

Tested versions: R 4.5.1; tidyverse 2.0; data.table 1.17; DESeq2 1.50;
arrow 22.

## Repository layout

```text
T21-eQTL/
  README.md                  # this file (canonical doc)
  LICENSE                    # BSD-2-Clause Plus Patent
  install_packages.R         # one-time CRAN + Bioconductor install
  environment.yml            # conda alternative
  download_gtex.sh           # fetch the GTEx v10 chr21 whole-blood allpairs parquet
  .github/workflows/ci.yml   # CI: parse-check scripts, run tests, upload coverage
  codecov.yml                # Codecov settings for the CI coverage upload
  .coderabbit.yaml           # CodeRabbit review settings (used on the fork)
  scripts/                   # production pipeline (see Pipeline section)
    lib/                     # shared helpers, unit-tested under tests/
    archive/                 # legacy + supplementary scripts (see docs/decisions.md)
  tests/
    testthat.R               # test runner: Rscript tests/testthat.R
    testthat/                # unit tests for scripts/lib/
  data/                      # inputs - mostly .gitignored
    HTP_WholeBlood_RNAseq_Counts_Synapse.txt   # 3.9 GB raw counts
    P4C_metadata_021921_Costello.txt           # sample metadata
    P4C_Comorbidity_020921.tsv                 # optional comorbidities
    chr21_ds_PASS.csv                          # 6.0 GB T21 chr21 genotypes
    chr21_ctrl_PASS.csv                        # 103 MB Control chr21 genotypes
    GTEx_Analysis_v10_..._Whole_Blood.v10.allpairs.chr21.parquet
                                               # GTEx allpairs (chr21)
    gtex/                                      # GTEx v10 median TPM per tissue (script 19 downloads it)
    raw/                                       # placeholder for raw downloads - .gitignored
    processed/                                 # script outputs - .gitignored
  config/runs/               # run definitions: baseline.R, adjusted.R
  results/
    runs/<name>/             # per-run outputs (.gitignored): processed/, tables/, figures/
    archive/2026-09-04_flat/ # pre-refactor flat outputs + drift_audit.md (historical)
    archive/                 # outputs of earlier pipeline versions (.gitignored)
  docs/
    summary.Rmd              # results summary; reads results/runs/<run>/ and renders to summary.md
    summary.md               # rendered summary (headline numbers and figures)
    decisions.md             # decision log, legacy notes, gotchas
    figures/                 # PNG copies of the pipeline figures embedded by summary.md
    package_installation_info.txt
```

## Pipeline

The production pipeline is the 00-07 chain plus 11 (eQTL figures), 15-18 (the eQTL test's power, base rate, negative control and effect-size comparison), 19 (GTEx tissue expression), 10 (run comparison) and the S1 supplement, run per named run (see Runs), in `scripts/`, backed by shared
helpers in `scripts/lib/` (`cohort.R` - analysis-cohort definition;
`chr21_threshold.R` - chr21-internal robust outlier test, annotation only;
`lane_rules.R` -
the sig_lane classification rule (`assign_sig_lane`); `eqtl_fit.R` -
vectorized per-variant regressions, the gene-level permutation test and
the runner scripts 03 and 16 share (`gene_level_tests`);
`eqtl_figures.R` - table preparation for script 11 (dosage-panel variants, best-variant effect sizes, chr21 map bands);
`allelic_fc.R` - allelic fold change fitted at ploidy 3, for scripts 15 and 18;
`inflation.R` - the distal-variant scan and genomic inflation, for script 17;
`gtex_expression.R` - the GTEx v10 median-TPM table (download, reader), the per-gene tissue summary, and counts to TPM over GENCODE exon-union lengths, for script 19;
`table1.R` - cohort-characteristics table helpers; `sankey_flow.R` - lane flow and its SankeyMATIC serialisation; `ploidy_distributions.R` - uncorrected vs ploidy-corrected log2FC tables for script 06; `biotypes.R` - the chr21 target biotype set, `TARGET_BIOTYPES`; `covariates.R` - CyTOF table reader, composition fractions, design formula, the expression-artifact adjustment and the per-covariate attribution; `run.R` - run configuration and paths; `karyotype_subtype.R` - INCLUDE subtype join).

| Script | Purpose | Output |
|---|---|---|
| 00_preprocess_data | Long -> wide gene x sample matrix; match metadata | `data/processed/count_matrix.csv`, `sample_metadata.csv`, `gene_annotations.csv` |
| 01_deseq2_analysis | Trisomy-aware DESeq2 (ploidy normalization matrix; chr21 excluded from size factors; betaPrior=FALSE) | `results/runs/<name>/tables/deseq2_all_genes_ploidy_normalized.csv`, `deseq2_chr21_genes_both_analyses.csv`, `deseq2_all_genes_both_analyses.csv`, QC PDFs |
| 02_filter_genotypes | Select deviating chr21 genes (`norm_padj < 0.01` AND `abs(norm_log2FC) >= log2(4/3)`, the tier-2 cut) among `TARGET_BIOTYPES` with baseMean >= 30, with the chr21-internal outlier annotation (`dev_z`, `q_outlier`); pull their GTEx allpairs cis variants (`pval_nominal <= 1e-4`), stream the PASS genotypes at those positions, and attach the minor-allele reference (GTEx `af`, gnomAD v4.1 AF) | `results/runs/<name>/processed/eqtl_supported_genes.csv`, `eqtl_target_variants.csv`, `genotypes_filtered.csv` |
| 03_t21_dosage_boxplots | Per-(variant, gene) within-T21 regressions of the run's expression artifact on genotype dosage, slopes also per minor allele; gene-level cis-eQTL permutation test (`scripts/lib/eqtl_fit.R`) | `results/runs/<name>/tables/t21_dosage_per_variant.csv`, `eqtl_allele_alignment.csv`, `eqtl_gene_level_perm.csv`, `t21_representative_variants.csv` |
| **04_chr21_lane_assignment** | **Per-gene lane assignment.** Hunter's padj rule with the two-tier magnitude cut (`norm_padj < 0.01` AND `abs(norm_log2FC) >= log2(4/3)`; tier 1 at log2(1.5)) is the classification split (`sig_lane`, applied by `scripts/lib/lane_rules.R`); deviating genes get an `eqtl_lane` terminal from the gene-level permutation test | `results/runs/<name>/tables/chr21_lane_assignments.csv`, `chr21_lane_summary.csv`, `chr21_k_sensitivity.csv` |
| 05_sankeymatic_export | Lane flow (Classification -> Sub-category -> eQTL terminal) serialised as SankeyMATIC input (`scripts/lib/sankey_flow.R`); the figure itself is rendered at https://sankeymatic.com/build/ | `results/runs/<name>/tables/chr21_lane_sankeymatic_input.txt`, `chr21_lane_flow.csv` |
| 06_chr21_distribution_panel | Density + ECDF of uncorrected vs ploidy-corrected log2FC for chr21 genes of the target biotypes, with chr22 on both scales as the control that the correction leaves unchanged (`scripts/lib/ploidy_distributions.R`) | `results/runs/<name>/figures/ploidy_correction_distributions.{pdf,png}`, `results/runs/<name>/tables/ploidy_correction_distribution_stats.csv` |
| 07_three_panel_figure | Two volcano figures: `volcano_all_genes` (A/B uncorrected and ploidy-corrected, all target-biotype genes with chr21 highlighted, no labels) and `volcano_chr21` (chr21 only, deviating genes labelled and coloured by direction, tier 1 bold); labels read from the lane table | `results/runs/<name>/figures/volcano_all_genes.{pdf,png}`, `volcano_chr21.{pdf,png}` |
| 11_eqtl_figures | eQTL-stage figures from existing tables: `eqtl_dosage_panels` (expression by dosage in T21 on each deviating gene's best variant, plotted on the allele whose GTEx effect runs the way the gene deviates; DE high and DE low panels), `chr21_deviating_map` (a band per deviating gene along chr21, coloured by eQTL outcome); writes the within-T21 slope per minor allele at each best variant, checked against the stored test p. Helpers in `scripts/lib/eqtl_figures.R` | `results/runs/<name>/figures/eqtl_dosage_panels.{pdf,png}`, `chr21_deviating_map.{pdf,png}`, `tables/eqtl_best_variant_effects.csv` |
| 15_spike_in_power | Power of the gene-level test for each tested deviating gene: its own T21 genotypes and cis variant set, its expression permuted across subjects, plus an allelic effect of a given aFC at a causal variant (ploidy-3 model), scored as script 03 scores it. `T21_N_SIM` (200) simulations per point | `results/runs/<name>/tables/spike_in_power.csv`, `spike_in_power_genes.csv`, `figures/spike_in_power.{pdf,png}` |
| 16_base_rate_expected_dosage | The identical gene-level test (same runner, variant universe, subjects, expression artifact and thresholds) on every Expected-dosage gene with GTEx cis variants: the reference set for the deviating genes | `results/runs/<name>/tables/base_rate_expected_dosage.csv`, `base_rate_summary.csv`, `processed/base_rate_genotypes.csv` |
| 17_distal_inflation | Negative control of the test: the run's T21 expression of every assessable chr21 gene (and chr22 genes as a reference) regressed on every common chr21 variant at least `distal_min_distance` (5 Mb) from the gene's TSS, with the run's covariates projected out of both sides; genomic inflation lambda and QQ against joint subject permutations (`scripts/lib/inflation.R`) | `results/runs/<name>/tables/distal_inflation_*.csv`, `figures/distal_inflation.{pdf,png}` |
| 18_detection_vs_effect_size | Deviating genes against the Expected-dosage GTEx eGenes (script 16) at matched effect size: detection expected from a logistic fit on GTEx \|aFC\| (Poisson-binomial) and the odds ratio; and each dataset's own top cis-eQTL effect, T21 aFC at the T21 best variant over GTEx aFC at its lead variant, deviating genes against the Expected-dosage eGenes in the same GTEx \|aFC\| range. Figure: that ratio (A) and spike-in power per deviating gene (B, script 15) | `results/runs/<name>/tables/detection_vs_effect_size{,_summary}.csv`, `top_variant_effect_ratio{,_summary}.csv`, `figures/eqtl_detection_overview.{pdf,png}` |
| 19_gtex_tissue_expression | Manuscript issue #30: is each chr21 gene meaningfully expressed in whole blood? Every classified chr21 gene against the GTEx v10 median TPM of 68 tissues (whole-blood TPM, rank of whole blood, top tissue, tissue specificity tau, percentile among all chr21 and among Expected-dosage genes), and against the eQTLs: GTEx whole-blood tested/eGene status and the within-T21 test (scripts 03, 16), whole-blood TPM compared between outcomes (Wilcoxon) a logistic curve of T21 detection on whole-blood TPM over the Expected-dosage genes, and the GTEx whole-blood eGene rate of the deviating genes against the other GTEx-tested chr21 genes (Fisher); the run's T21 participants on the TPM scale (counts over GENCODE v44 exon-union lengths) against GTEx whole blood for every chr21 gene, with the genome-wide HTP/GTEx offset (HTP is globin-depleted). Figure: A whole-blood TPM by lane, B GTEx eQTL strength against TPM for the chr21 genes GTEx tested, shaped by the T21 test outcome (as in D), C whole blood among the 68 tissues (coloured by organ-system group) per deviating gene, D the T21 participants' median TPM against GTEx for every chr21 gene. Annotation only | `results/runs/<name>/tables/gtex_tissue_expression_{genes,summary,eqtl,detection,scale,egene_rate}.csv`, `figures/gtex_tissue_expression.{pdf,png}` |
| 10_compare_runs | Per-gene comparison of two runs' lane tables (lane, tier, norm_log2FC, eQTL call, attenuation) and a lane-transition table; `--run adjusted --against baseline` | `results/runs/adjusted/tables/run_comparison_adjusted_vs_baseline.csv`, `lane_transitions_adjusted_vs_baseline.csv`, `figures/run_comparison_adjusted_vs_baseline.{pdf,png}` |
| S1_covariate_evidence | Supplement per run: A covariates vs karyotype; B cell-fraction reach into expression (composition runs); C mosaic subtype vs chr21 index and genome-wide expression; D per-covariate attribution of the baseline-to-adjusted fold-change shift with bootstrap intervals | `results/runs/<name>/tables/S1_*.csv`, `figures/S1_covariate_evidence.{pdf,png}` |
| audit_baseline_drift | Table-by-table comparison of a regenerated run against the archived 2026-09-04 flat outputs; strict on DESeq2 tables and lane classification, reports expression-scale-dependent tables | `results/archive/2026-09-04_flat/drift_audit.md` |

Legacy and supplementary scripts live under `scripts/archive/`; the
catalog is in [docs/decisions.md](docs/decisions.md).


## Outputs

### Tables

`results/runs/<name>/tables/`:
- `deseq2_chr21_genes_both_analyses.csv` - chr21 DESeq2 results (raw + ploidy-
  corrected in one row per gene).
- `deseq2_all_genes_ploidy_normalized.csv` - genome-wide ploidy-corrected
  DESeq2 results (reference; classification uses the chr21 table directly).
- **`chr21_lane_assignments.csv`** - canonical per-gene lane table (read
  this for the headline numbers).
- `chr21_lane_summary.csv` - lane counts (all chr21 + after paper filters).
- `chr21_lane_flow.csv` - lane-flow paths (level2 / level3 / level4 counts) behind the SankeyMATIC input.
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

`results/runs/<name>/figures/`:
- Lane-flow Sankey: rendered at https://sankeymatic.com/build/ from
  `results/runs/<name>/tables/chr21_lane_sankeymatic_input.txt` (script 05); the
  tracked render is `docs/figures/Sankey.png`.
- `ploidy_correction_distributions.{pdf,png}` - uncorrected vs
  ploidy-corrected log2FC distributions for chr21, with chr22 as the
  unchanged control (script 06).
- `volcano_all_genes.{pdf,png}`, `volcano_chr21.{pdf,png}` - volcano figures (script 07).
- `eqtl_dosage_panels.{pdf,png}`, `chr21_deviating_map.{pdf,png}` -
  eQTL-stage figures (script 11).
- `eqtl_detection_overview.{pdf,png}` - top cis-eQTL effect in T21 against
  GTEx for deviating and matched Expected-dosage genes, and spike-in power
  per deviating gene (scripts 15, 16, 18).
- `distal_inflation.{pdf,png}` - the distal-variant negative control (script 17).
- `gtex_tissue_expression.{pdf,png}` - GTEx v10 whole-blood expression of the
  chr21 genes by lane, against GTEx whole-blood eQTL strength and against the
  T21 participants' TPM; per deviating gene, whole blood among the GTEx
  tissues (script 19).


## Further documentation

[docs/decisions.md](docs/decisions.md) - decision log (retired
classification filters, the replaced eQTL rule, the retired
neighborhood check), legacy data sources and terminology, the `scripts/archive/`
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
