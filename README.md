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
Rscript scripts/11_eqtl_figures.R --run baseline                # eQTL dosage panels, effect sizes, chr21 map
Rscript scripts/13_af_deviation_bound.R --run baseline         # can the cis-eQTLs move the group ratio? (AF bound)
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

- `baseline`: Hunter-style, design `~ karyotype`, no exclusions. The primary
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
  config/runs/               # run definitions: baseline.R, adjusted.R
  results/
    runs/<name>/             # per-run outputs (.gitignored): processed/, tables/, figures/
    archive/2026-09-04_flat/ # pre-refactor flat outputs + drift_audit.md (historical)
    archive/                 # outputs of earlier pipeline versions (.gitignored)
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
positions of the eQTL universe. Derived from the HTP VCFs on the INCLUDE
Data Hub; `scripts/drs_download.R` and `scripts/download_vcf_array.sh`
fetch those VCFs from a DRS manifest. The VCF-to-PASS-CSV step is not in
this repository.

### 5. GTEx whole-blood eQTLs

Chr21 extract from GTEx v10 allpairs (every cis-window variant tested per
gene, regardless of significance). Filename:
`data/GTEx_Analysis_v10_QTLs_GTEx_Analysis_v10_eQTL_all_associations_Whole_Blood.v10.allpairs.chr21.parquet`.
Script 02 applies `pval_nominal <= 1e-4` to keep the variant universe
manageable: a fixed stand-in for GTEx's per-gene nominal threshold (chr21
whole-blood median 1.4e-4) that reproduces the signif_pairs eGene set. The
all-tissue GTEx files under `data/` are not read by any script.

GTEx distributes the v10 all-associations results per tissue and per
chromosome as parquet files in a requester-pays Google Cloud bucket
(`gs://gtex-resources/GTEx_Analysis_v10_QTLs/GTEx_Analysis_v10_eQTL_all_associations/`;
see https://gtexportal.org/home/downloads/adult-gtex/qtl). The chr21
whole-blood file is used exactly as distributed: no conversion step, and
no python or pyarrow, is involved. To re-download it, run from the repo
root

```bash
GCP_BILLING_PROJECT=<your-gcp-project> bash download_gtex.sh
```

which fetches `Whole_Blood.v10.allpairs.chr21.parquet` with `gcloud
storage cp` (or `gsutil cp`), billing the egress to the named project,
validates it with `arrow` (a parquet carrying the columns script 02 reads:
`gene_id`, `variant_id`, `pval_nominal`), and only then moves it to the
filename above. Without the Google Cloud SDK the script prints the object
path for a manual download and exits non-zero; without R and `arrow` it
keeps the unvalidated download as a `.part` file and exits non-zero.

### 6. HTP cell-type data (Synapse CyTOF; read by script 01 in composition runs)

**Content**: HTP mass cytometry (CyTOF) of CD45+ CD66low white blood cells,
FlowSOM-clustered into 20 cell populations, as percent of CD45+ CD66low
cells per sample (`HTP_CyTOF_CD45posCD66low_FlowSOM_cluster_percentage_Synapse.txt`,
long format: `LabID`, `Cell_cluster_name`, `Value`), with a column
dictionary and the antibody panel. 388 samples, one per subject. 372 of the
397 analysis-cohort samples have an exact `LabID` match (282 T21, 90
Control). The CD66low gate excludes granulocytes, so neutrophil fraction,
the dominant axis of whole-blood RNA-seq composition, is not measured here.

The raw single-cell data behind that table is FlowRepository experiment
FR-FCM-Z5GE ("Human Trisome Project CyTOF CD45pos CD66low", Araya /
Galbraith, PMID 37379383): 388 FCS files, one per sample, each a
batch-adjusted 10,000-event subset of the same CD45+ CD66low gate, plus the
same antibody panel. Kept as `data/flowrepository/FlowRepository_FR-FCM-Z5GE_files.zip`
(704 MB, gitignored) for provenance; the Synapse percentage table is its
processed output and is what any analysis should read. It adds no
granulocyte information.

Synapse folder `syn31488783`, downloaded to `data/synapse/syn31488783/` by

```bash
printf '%s' '<synapse personal access token>' > ~/.synapse_token && chmod 600 ~/.synapse_token
bash download_synapse_celltypes.sh
```

which runs the Synapse Python client from a virtualenv it creates under
`~/.virtualenvs/synapser`. The token is read from `SYNAPSE_AUTH_TOKEN` or
`~/.synapse_token` and never stored in the repo. Read by script 01 in composition runs
(`config/runs/adjusted.R`) as 19 cell-fraction covariates.

### 7. Karyotype subtype (INCLUDE Data Hub; read by script 00)

**Files**: `data/karyotype/simple_which_T21.csv` (INCLUDE participant id ->
`genotype`: `T21`, `DS_T21`, `mosaic_T21`, `translocation_T21`, `D21`; 622
participants) and `data/karyotype/personnameswitch.csv` (INCLUDE participant
id -> HTP external participant id; 835 rows, 578 with an HTP id). Joined by
subject to the analysis cohort in `data/processed/karyotype_subtype.csv`
(script 00 writes it and adds `karyotype_subtype` to `sample_metadata.csv`). All 302
T21 subjects map: 252 `T21`, 32 `DS_T21`, 9 `mosaic_T21`, 9
`translocation_T21`; all 95 controls are `D21`. The per-sample chr21
expression index (`results/runs/<name>/tables/chr21_dosage_index_per_sample.csv`,
median ratio to the control median over well-expressed chr21 genes) agrees:
mosaic median 1.07 against 1.35 for full trisomy. The baseline run keeps every T21
subject at ploidy 1.5, mosaics included; the adjusted run excludes the mosaics.

## Pipeline

The production pipeline is the 00-07 chain plus 11 (eQTL figures), 10 (run comparison) and the S1 supplement, run per named run (see Runs), in `scripts/`, backed by shared
helpers in `scripts/lib/` (`cohort.R` - analysis-cohort definition;
`chr21_threshold.R` - chr21-internal robust outlier test, annotation only;
`lane_rules.R` -
the sig_lane classification rule (`assign_sig_lane`); `eqtl_fit.R` -
vectorized per-variant regressions and the gene-level permutation test;
`eqtl_controls.R` - the shared gene-level test runner and the standalone
negative (decoy variant set) and positive (GTEx eGene) controls;
`eqtl_figures.R` - table preparation for script 11 (dosage-panel variants, per-gene rows, best-variant effect sizes, chr21 map bands);
`table1.R` - cohort-characteristics table helpers; `sankey_flow.R` - lane flow and its SankeyMATIC serialisation; `ploidy_distributions.R` - uncorrected vs ploidy-corrected log2FC tables for script 06; `biotypes.R` - the chr21 target biotype set, `TARGET_BIOTYPES`; `covariates.R` - CyTOF table reader, composition fractions, design formula, the expression-artifact adjustment and the per-covariate attribution; `run.R` - run configuration and paths; `karyotype_subtype.R` - INCLUDE subtype join).

| Script | Purpose | Output |
|---|---|---|
| 00_preprocess_data | Long -> wide gene x sample matrix; match metadata | `data/processed/count_matrix.csv`, `sample_metadata.csv`, `gene_annotations.csv` |
| 01_deseq2_analysis | Trisomy-aware DESeq2 (ploidy normalization matrix; chr21 excluded from size factors; betaPrior=FALSE) | `results/runs/<name>/tables/deseq2_all_genes_ploidy_normalized.csv`, `deseq2_chr21_genes_both_analyses.csv`, `deseq2_all_genes_both_analyses.csv`, QC PDFs |
| 02_filter_genotypes | Select deviating chr21 genes by Hunter et al.'s rule (`norm_padj < ALPHA_DE` AND `abs(norm_log2FC) >= DEVIATION_LFC`), restrict to `TARGET_BIOTYPES` (protein-coding, lncRNA, pseudogene: the biotypes GTEx tests; `scripts/lib/biotypes.R`), apply Hunter's baseMean >= 30 coverage floor, compute the chr21-internal outlier annotation (`dev_z`, `q_outlier`), pull GTEx allpairs cis variants, stream PASS files for those positions; add the `n_positive_controls` strongest GTEx whole-blood eGenes among expressed, non-repeat, non-deviating chr21 genes to the roster as `gene_set == "positive_control"` (`scripts/lib/eqtl_controls.R`) | `data/processed/eqtl_supported_genes.csv`, `eqtl_target_variants.csv`, `genotypes_filtered.csv` |
| 03_t21_dosage_boxplots | Per-(variant, gene) within-T21 regressions of the run's expression artifact (log2-CPM, covariate effects removed) on alt-allele dosage; gene-level cis-eQTL permutation test (`scripts/lib/eqtl_fit.R`); standalone controls of that test: negative (each deviating gene against the cis variants of a distant tested gene) and positive (the GTEx eGenes from script 02), run through the same runner (`scripts/lib/eqtl_controls.R`) | `results/runs/<name>/tables/t21_dosage_per_variant.csv`, `t21_representative_variants.csv`, `eqtl_gene_level_perm.csv`, `eqtl_control_negative.csv`, `eqtl_control_positive.csv`, `eqtl_controls_summary.csv` |
| **04_chr21_lane_assignment** | **Per-gene lane assignment.** Hunter's padj rule with the two-tier magnitude cut (`norm_padj < 0.01` AND `abs(norm_log2FC) >= log2(4/3)`; tier 1 at log2(1.5)) is the classification split (`sig_lane`, applied by `scripts/lib/lane_rules.R`); deviating genes get an `eqtl_lane` terminal from the gene-level permutation test | `results/runs/<name>/tables/chr21_lane_assignments.csv`, `chr21_lane_summary.csv`, `chr21_k_sensitivity.csv` |
| 05_sankeymatic_export | Lane flow (Classification -> Sub-category -> eQTL terminal) serialised as SankeyMATIC input (`scripts/lib/sankey_flow.R`); the figure itself is rendered at https://sankeymatic.com/build/ | `results/runs/<name>/tables/chr21_lane_sankeymatic_input.txt`, `chr21_lane_flow.csv` |
| 06_chr21_distribution_panel | Density + ECDF of uncorrected vs ploidy-corrected log2FC for chr21 genes of the target biotypes, with chr22 on both scales as the control that the correction leaves unchanged (`scripts/lib/ploidy_distributions.R`) | `results/runs/<name>/figures/ploidy_correction_distributions.{pdf,png}`, `results/runs/<name>/tables/ploidy_correction_distribution_stats.csv` |
| 07_three_panel_figure | Two volcano figures: `volcano_all_genes` (A/B uncorrected and ploidy-corrected, all target-biotype genes with chr21 highlighted, no labels) and `volcano_chr21` (chr21 only, deviating genes labelled and coloured by direction, tier 1 bold); labels read from the lane table | `results/runs/<name>/figures/volcano_all_genes.{pdf,png}`, `volcano_chr21.{pdf,png}` |
| 11_eqtl_figures | eQTL-stage figures from existing tables: `eqtl_dosage_panels` (expression by alt dosage in T21, one panel per gene on its best variant, in DE high / DE low / positive-control / negative-control groups), `eqtl_effect_sizes` (scatter of the within-T21 slope at the best variant against -log10 gene-level permutation q, for deviating genes by direction, their decoy variant sets and the positive controls, with the q = 0.05 line; refit p checked against the stored test p; writes `tables/eqtl_best_variant_effects.csv`), `chr21_deviating_map` (schematic chr21 with a band per deviating gene coloured by eQTL outcome, labels by direction; positions from the roster TSS and `data/chr21_gene_positions.csv`). Helpers in `scripts/lib/eqtl_figures.R` | `results/runs/<name>/figures/eqtl_dosage_panels.{pdf,png}`, `eqtl_effect_sizes.{pdf,png}`, `chr21_deviating_map.{pdf,png}`, `tables/eqtl_best_variant_effects.csv` |
| 10_compare_runs | Per-gene comparison of two runs' lane tables (lane, tier, norm_log2FC, eQTL call, attenuation) and a lane-transition table; `--run adjusted --against baseline` | `results/runs/adjusted/tables/run_comparison_adjusted_vs_baseline.csv`, `lane_transitions_adjusted_vs_baseline.csv`, `figures/run_comparison_adjusted_vs_baseline.{pdf,png}` |
| S1_covariate_evidence | Supplement per run: A covariates vs karyotype; B cell-fraction reach into expression (composition runs); C mosaic subtype vs chr21 index and genome-wide expression; D per-covariate attribution of the baseline-to-adjusted fold-change shift with bootstrap intervals | `results/runs/<name>/tables/S1_*.csv`, `figures/S1_covariate_evidence.{pdf,png}` |
| audit_baseline_drift | Table-by-table comparison of a regenerated run against the archived 2026-09-04 flat outputs; strict on DESeq2 tables and lane classification, reports expression-scale-dependent tables | `results/archive/2026-09-04_flat/drift_audit.md` |

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
- Flag low-expression genes with Hunter et al.'s minimum coverage floor
  (baseMean < 30), applied identically in scripts 02 and 04.

After ploidy normalization, the appropriate null on chr21 is back to `FC = 1`,
so DESeq2's standard p-value testing applies cleanly.



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
- `eqtl_dosage_panels.{pdf,png}`, `eqtl_effect_sizes.{pdf,png}`,
  `chr21_deviating_map.{pdf,png}` - eQTL-stage figures with the controls (script 11).


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
