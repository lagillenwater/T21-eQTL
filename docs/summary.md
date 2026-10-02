T21-eQTL results summary
================
2026-10-02

Rendered at commit `007f38d` (the commit checked out when this document
was last knitted; `results/` is gitignored, so the outputs themselves
are not versioned). This document reads:

- `results/runs/baseline/tables/`: `chr21_lane_assignments.csv` (script
  04\) - every number and gene list in the baseline sections;
  `eqtl_controls_summary.csv`, `eqtl_control_positive.csv` and
  `eqtl_control_negative.csv` (script 03);
  `positive_control_matching.csv` (script 02);
  `eqtl_best_variant_effects.csv` (script 11).
- `config/runs/baseline.R` - the control thresholds quoted in the text.
- `results/runs/adjusted/`: `processed/cohort_roster.csv` (script 01);
  `tables/chr21_lane_assignments.csv` and `chr21_lane_summary.csv`
  (script 04); `tables/eqtl_controls_summary.csv` (script 03);
  `tables/lane_transitions_adjusted_vs_baseline.csv` and
  `run_comparison_adjusted_vs_baseline.csv` (script 10);
  `tables/S1_attribution.csv` and `S1_covariates_vs_karyotype.csv` (S1).
- `docs/figures/` - the PNGs embedded under Figures and in the
  adjusted-run section (copies of the outputs of scripts 06, 07, 10, 11
  and S1), and `Sankey.png`, the SankeyMATIC render of
  `results/runs/<run>/tables/chr21_lane_sankeymatic_input.txt` (script
  05).

To refresh after a pipeline run:
`Rscript -e 'rmarkdown::render("docs/summary.Rmd")'`. Methodology and
its history: `README.md` and `docs/decisions.md`.

## Overview

Of the **318** chr21 genes classified (**160** protein-coding, **129**
lncRNA, **29** pseudogene: the biotypes GTEx whole blood tests for
cis-eQTLs): **130** are Expected dosage, **12** are flagged as high
repeats and **153** low-expression (baseMean below Hunter et al.’s
coverage floor of 30), and **23** deviate (**10** DE_high, **13**
DE_low; **7** tier 1 (padj \< 0.01 and \|corrected log2FC\| \>=
log2(1.5), a fold-change magnitude of 1.5 in either direction), **16**
tier 2, near-threshold (padj \< 0.01 and log2(4/3) \<= \|corrected
log2FC\| \< log2(1.5))). Of the deviating genes, **14** have a
detectable GTEx whole-blood cis-eQTL, **6** were tested and have none,
and **3** have no GTEx cis coverage.

| sig_lane        | cis_eqtl | no_GTEx_data | no_cis_eqtl | not_evaluated | Total |
|:----------------|---------:|-------------:|------------:|--------------:|------:|
| Low_expression  |        0 |            0 |           0 |           153 |   153 |
| Expected_dosage |        0 |            0 |           0 |           130 |   130 |
| DE_low          |        6 |            1 |           6 |             0 |    13 |
| High_repeats    |        0 |            0 |           0 |            12 |    12 |
| DE_high         |        8 |            2 |           0 |             0 |    10 |

Lane counts by eQTL terminal. cis_eqtl is a detection result (q_gene_bh
\< 0.05), not an ‘explained by eQTL’ claim.

## Deviating genes

| Gene         | lane    | tier | log2FC | padj    | eqtl         | q_perm |
|:-------------|:--------|-----:|-------:|:--------|:-------------|-------:|
| ABCC13       | DE_high |    1 |   1.45 | 2.2e-22 | cis_eqtl     | 0.0220 |
| TSPEAR       | DE_high |    1 |   0.75 | 1.3e-09 | cis_eqtl     | 0.0029 |
| RIPK4        | DE_high |    1 |   0.64 | 7.3e-09 | cis_eqtl     | 0.0067 |
| AP001610.2   | DE_high |    1 |   0.64 | 3.1e-05 | no_GTEx_data |     NA |
| CYYR1        | DE_high |    2 |   0.57 | 3.6e-09 | cis_eqtl     | 0.0410 |
| COL6A2       | DE_high |    2 |   0.56 | 4.2e-12 | cis_eqtl     | 0.0150 |
| YBEY         | DE_high |    2 |   0.50 | 8.7e-24 | cis_eqtl     | 0.0029 |
| ERG          | DE_high |    2 |   0.50 | 4.9e-03 | no_GTEx_data |     NA |
| MX1          | DE_high |    2 |   0.47 | 1.7e-04 | cis_eqtl     | 0.0140 |
| AATBC        | DE_high |    2 |   0.44 | 2.1e-14 | cis_eqtl     | 0.0220 |
| OLIG2        | DE_low  |    1 |  -1.19 | 5.0e-19 | no_cis_eqtl  | 0.1700 |
| COL6A1       | DE_low  |    1 |  -0.83 | 2.3e-17 | cis_eqtl     | 0.0029 |
| TEKT4P2      | DE_low  |    1 |  -0.64 | 4.2e-04 | no_cis_eqtl  | 0.1500 |
| AF165147.1   | DE_low  |    2 |  -0.57 | 1.4e-11 | no_cis_eqtl  | 0.1500 |
| PDE9A        | DE_low  |    2 |  -0.55 | 3.4e-10 | cis_eqtl     | 0.0029 |
| ZBTB21       | DE_low  |    2 |  -0.51 | 1.7e-33 | no_GTEx_data |     NA |
| RBM11        | DE_low  |    2 |  -0.47 | 4.0e-06 | no_cis_eqtl  | 0.0590 |
| LINC01679    | DE_low  |    2 |  -0.45 | 4.8e-08 | cis_eqtl     | 0.0029 |
| PAXBP1-AS1   | DE_low  |    2 |  -0.45 | 8.3e-17 | no_cis_eqtl  | 0.8100 |
| PCBP3        | DE_low  |    2 |  -0.43 | 7.0e-09 | cis_eqtl     | 0.0029 |
| BACE2        | DE_low  |    2 |  -0.43 | 1.1e-09 | no_cis_eqtl  | 0.1500 |
| C21orf62-AS1 | DE_low  |    2 |  -0.42 | 1.2e-10 | cis_eqtl     | 0.0050 |
| ICOSLG       | DE_low  |    2 |  -0.42 | 4.4e-03 | cis_eqtl     | 0.0029 |

Ploidy-corrected stats and the gene-level permutation q. eqtl is a
detection result; it does not identify what produces the deviation.

Deviating genes that carry no detected cis-eQTL. These deviations are
not accompanied by detectable cis-regulatory variation; the pipeline
does not attribute them to a mechanism:

- Tested against GTEx, no detection (6): **AF165147.1, BACE2, OLIG2,
  PAXBP1-AS1, RBM11, TEKT4P2**
- No GTEx cis coverage, untestable here (3): **AP001610.2, ERG, ZBTB21**

## Controls for the eQTL test

Two standalone gene sets go through the identical gene-level permutation
test (same runner, 1000 permutations, BH within the set, q \< 0.05).

**Negative control, unlinked variants.** Each deviating gene’s
expression is tested against 5 cis variant sets of other tested genes
(deviating genes and positive controls) whose TSS is at least 5 Mb away,
so no decoy set can be in LD with the gene’s own locus. Sets with the
closest variant count are taken first, and decoy variants below GTEx MAF
0.05 are dropped. Same genotypes, same subjects. Every decoy test is
null, so BH should detect none (with no true effects it holds the chance
of any detection to the FDR level), and about 5% of tests should reach
nominal permutation p \< 0.05.

**Positive control, matched GTEx eGenes.** Up to one control per tested
deviating gene, matched in order of decreasing GTEx effect size until
the one-per-locus rule below has used up the candidates (here 13 of the
20 tested genes), drawn from the GTEx whole-blood eGenes (q \< 0.0001,
with at least 10 variants at the GTEx pval cut) among expressed,
non-repeat, non-deviating chr21 genes: the candidate nearest the
deviating gene in GTEx allelic fold change and expression level, one per
locus (TSS at least 1 Mb from every other control and 100 kb from any
deviating gene), never a gene in a GRCh38 false duplication. Each is
tested on its own cis variants. Most should be detected; if not, the
test lacks power for effects of the size it is asked to find.

| run | set | n_tested | n_detected | pct_detected | n_nominal_p05 | expectation |
|:---|:---|---:|---:|---:|---:|:---|
| baseline | observed_deviating | 20 | 14 | 70.0 | 15 | the result |
| baseline | negative_unlinked_variants | 100 | 0 | 0.0 | 10 | about 0 at q; about 5% at nominal p \< 0.05 |
| baseline | positive_gtex_egenes | 13 | 11 | 84.6 | 11 | most detected |
| adjusted | observed_deviating | 10 | 7 | 70.0 | 8 | the result |
| adjusted | negative_unlinked_variants | 50 | 0 | 0.0 | 6 | about 0 at q; about 5% at nominal p \< 0.05 |
| adjusted | positive_gtex_egenes | 10 | 8 | 80.0 | 8 | most detected |

eQTL test: observed set and the two standalone controls, both runs.

| Control | Matched to | GTEx abs aFC | matched gene abs aFC | GTEx min p | n variants | q | detected |
|:---|:---|---:|---:|:---|---:|---:|:---|
| CYP4F29P | CYYR1 | 2.76 | 3.07 | 5.6e-50 | 291 | 0.0014 | TRUE |
| LINC00189 | TEKT4P2 | 1.53 | 1.46 | 3.7e-21 | 610 | 0.0014 | TRUE |
| LRRC3 | ABCC13 | 1.64 | 1.39 | 4.0e-33 | 124 | 0.0014 | TRUE |
| TMPRSS3 | LINC01679 | 0.78 | 1.30 | 8.0e-20 | 107 | 0.0014 | TRUE |
| CBR3 | ICOSLG | 1.04 | 1.08 | 6.2e-15 | 201 | 0.0014 | TRUE |
| EVA1C | OLIG2 | 0.66 | 1.04 | 1.5e-22 | 102 | 0.1200 | FALSE |
| GET1 | YBEY | 0.57 | 0.93 | 5.0e-68 | 661 | 0.0014 | TRUE |
| KCNJ15 | MX1 | 0.56 | 0.82 | 1.4e-37 | 545 | 0.0014 | TRUE |
| ITSN1 | PDE9A | 0.50 | 0.72 | 8.4e-31 | 196 | 0.0014 | TRUE |
| DIP2A | COL6A2 | 0.71 | 0.65 | 1.3e-114 | 1005 | 0.0014 | TRUE |
| ADAMTS1 | TSPEAR | 0.37 | 0.63 | 1.9e-08 | 23 | 0.0026 | TRUE |
| MRPL39 | BACE2 | 0.29 | 0.60 | 7.3e-21 | 301 | 0.0140 | TRUE |
| NRIP1 | RIPK4 | 0.25 | 0.58 | 2.9e-08 | 14 | 0.3800 | FALSE |

Baseline positive control: one matched GTEx whole-blood eGene per tested
deviating gene, tested within T21 on its own cis variants. aFC = GTEx
allelic fold change (absolute log2).

| Gene         | decoy sets | detected (q \< 0.05) | nominal p \< 0.05 | smallest q |
|:-------------|-----------:|---------------------:|------------------:|-----------:|
| AATBC        |          5 |                    0 |                 0 |       0.60 |
| ABCC13       |          5 |                    0 |                 2 |       0.46 |
| AF165147.1   |          5 |                    0 |                 0 |       0.60 |
| BACE2        |          5 |                    0 |                 0 |       0.65 |
| C21orf62-AS1 |          5 |                    0 |                 0 |       0.91 |
| COL6A1       |          5 |                    0 |                 1 |       0.34 |
| COL6A2       |          5 |                    0 |                 2 |       0.32 |
| CYYR1        |          5 |                    0 |                 0 |       0.60 |
| ICOSLG       |          5 |                    0 |                 2 |       0.25 |
| LINC01679    |          5 |                    0 |                 0 |       0.65 |
| MX1          |          5 |                    0 |                 0 |       0.65 |
| OLIG2        |          5 |                    0 |                 0 |       0.60 |
| PAXBP1-AS1   |          5 |                    0 |                 0 |       0.65 |
| PCBP3        |          5 |                    0 |                 1 |       0.25 |
| PDE9A        |          5 |                    0 |                 0 |       0.60 |
| RBM11        |          5 |                    0 |                 0 |       0.60 |
| RIPK4        |          5 |                    0 |                 1 |       0.32 |
| TEKT4P2      |          5 |                    0 |                 0 |       0.60 |
| TSPEAR       |          5 |                    0 |                 1 |       0.27 |
| YBEY         |          5 |                    0 |                 0 |       0.80 |

Baseline negative control, per deviating gene: its expression tested
against each of its decoy variant sets. A gene with fewer sets had fewer
donors far enough away.

## Figures

Copies of the pipeline figures, tracked under `docs/figures/`. To
refresh them, re-run scripts 06, 07, 11 and S1 for each run (and
re-export the lane flow from SankeyMATIC using
`results/runs/<run>/tables/chr21_lane_sankeymatic_input.txt`), copy the
PNGs into `docs/figures/`, and re-render this document.

### Effect of ploidy correction on the log2FC distribution: chr21 uncorrected vs corrected, with chr22 as the unchanged control (script 06):

<img src="./figures/baseline_ploidy_correction_distributions.png" alt="Two-panel figure. Left: overlaid density curves of log2 fold change (T21 vs Control) for chr21 and chr22 genes (protein-coding, lncRNA, pseudogene), each on the uncorrected and the ploidy-corrected scale. The uncorrected chr21 curve is centred near log2(1.5); the ploidy-corrected chr21 curve is shifted left by that amount, centres on zero and overlaps chr22. The two chr22 curves coincide. Right: the corresponding empirical cumulative distributions. Vertical reference lines mark 0 and log2(1.5)." width="100%" />

### Volcano plots, all genes before and after ploidy correction (script 07):

<img src="./figures/baseline_volcano_all_genes.png" alt="Two volcano plots of log2 fold change (T21 vs Control) against minus log10 adjusted p-value, all genes of the target biotypes in grey with chr21 genes in purple. Left, uncorrected: the chr21 cloud sits to the right of zero around the 1.5-fold expectation, marked by a dotted line, and is almost uniformly significant. Right, ploidy-corrected: the chr21 cloud is centred on zero and most of its significance is gone." width="100%" />

### Volcano plots, chr21 only, deviating genes labelled (script 07):

<img src="./figures/baseline_volcano_chr21.png" alt="Two volcano plots restricted to chr21 genes, uncorrected on the left and ploidy-corrected on the right. Deviating genes are labelled and coloured by direction, red for higher than expected and blue for lower than expected, with tier 1 genes in bold; all other chr21 genes are grey." width="100%" />

### cis-eQTL effect size per gene, with its controls (script 11):

<img src="./figures/baseline_eqtl_effect_sizes.png" alt="Scatter plot with the within-T21 slope of expression per copy of the minor allele at the best variant on the x-axis and minus log10 of the gene-level permutation q on the y-axis, with a dashed line at q = 0.05. Green triangles are the matched positive-control GTEx eGenes, most on or near the top row with slopes from about -0.25 to 0.6, two below the line. Red and blue circles are deviating genes higher and lower than expected; most sit above the dashed line with slopes within about 0.45 of zero, and six sit below it. Grey diamonds are the decoy variant sets, all near the bottom, well below the line. The subtitle names the deviating genes with no GTEx variants." width="90%" />

Each point is one gene’s best variant: its within-T21 slope on the
x-axis (per copy of the minor allele) and the gene-level permutation q
on the y-axis, the value that decides detection. The dashed line is q =
0.05. The median absolute slope is 0.22 log2-CPM per allele on the
deviating genes’ own variants, 0.09 on their decoy sets, and 0.23 for
the positive controls. With 1000 permutations q cannot fall much below
0.001, so the strongest genes share the top row and differ only in
slope. Genes with no GTEx variants have neither value and are named in
the subtitle. Every slope is the best of many variants, so it is
inflated by that selection; the decoys show that a large best-variant
slope can occur with no cis link.

### Expression by genotype for every tested deviating gene (script 11; the controls are drawn the same way in `eqtl_dosage_controls`):

<img src="./figures/baseline_eqtl_dosage_panels.png" alt="Two rows of box-and-jitter panels, one per tested deviating gene, showing expression in T21 subjects against the dosage (0 to 3 copies) at the gene's best variant of the allele whose GTEx effect runs the way the gene deviates. Row A, DE high genes in red, is drawn on the expression-raising allele, so a panel that reproduces GTEx trends up; row B, DE low genes in blue, on the expression-lowering allele, so it trends down. Each strip names the gene, the permutation q, the variant, the plotted allele, the minor allele and its MAF, the GTEx direction and whether the within-T21 trend agrees. Every trend line is green, meaning the within-T21 trend matches the GTEx direction." width="100%" />

### Deviating genes along chr21 (script 11):

<img src="./figures/baseline_chr21_deviating_map.png" alt="Schematic of chromosome 21 as a horizontal bar with the centromere in dark grey. A coloured band marks the transcription start site of each deviating gene: green for cis-eQTL detected, orange for tested without detection, grey for no GTEx coverage. Gene names sit above the bar in red for genes higher than expected and below the bar in blue for genes lower than expected." width="100%" />

### Sankey plot generated with SankeyMATIC

<img src="./figures/Sankey.png" alt="Sankey diagram of the chr21 genes (protein-coding, lncRNA, pseudogene). From the full set, flows split into Expected dosage, Not assessable (which divides into High repeats and Low expression), and Outside dosage expectation. The last divides into DE high and DE low, which then terminate at cis-eQTL detected, no GTEx cis-eQTL data (labelled 'no GTEX QTL' in the figure), or no cis-eQTL detected. Node labels carry the gene counts." width="100%" />

## Covariate-adjusted run

The baseline result above follows Hunter et al.: no covariates. Age,
BMI, sample source and visit differ between T21 and control subjects in
this cohort (S1 panel A), blood cell composition differs by karyotype
and tracks thousands of genes (S1 panel B), and 9 T21 subjects are
mosaic, with a chr21 expression index near 1.07 rather than 1.35 (S1
panel C). The `adjusted` run (`config/runs/adjusted.R`) therefore fits
`~ age + sex + BMI + sample source + 19 CyTOF cell fractions + karyotype`
on the 354 samples (274 T21, 80 Control) that have every covariate, with
mosaic T21 excluded, and its eQTL stage reads the covariate-adjusted
expression artifact. It is reported beside the baseline, not instead of
it: BMI and composition are plausibly downstream of trisomy, so this run
estimates the deviation net of those pathways.

| reason                               |   n |
|:-------------------------------------|----:|
| exclude:karyotype_subtype=mosaic_T21 |   9 |
| keep_rule:no_wgs                     |   2 |
| missing:BMI                          |  12 |
| missing:cytof                        |  22 |

Samples excluded from the adjusted run, by reason.

| sig_lane             | eqtl_lane     | n_genes |
|:---------------------|:--------------|--------:|
| DE_high              | cis_eqtl      |       4 |
| DE_high              | no_GTEx_data  |       1 |
| DE_high              | no_cis_eqtl   |       1 |
| DE_low               | cis_eqtl      |       3 |
| DE_low               | no_GTEx_data  |       2 |
| DE_low               | no_cis_eqtl   |       2 |
| Expected_dosage      | not_evaluated |     129 |
| High_repeats         | not_evaluated |      12 |
| Low_expression       | not_evaluated |     153 |
| Not_DE_outside_noise | not_evaluated |      11 |

Adjusted run: chr21 genes by classification lane and eQTL terminal.
Not_DE_outside_noise = clears the magnitude cut but not padj \< 0.01 in
the adjusted model.

Of the 23 genes deviating in the baseline, 6 still deviate after
adjustment (OLIG2, PCBP3, BACE2, YBEY, COL6A2, ABCC13), 17 do not
(COL6A1, TEKT4P2, AF165147.1, PDE9A, ZBTB21, RBM11, LINC01679,
PAXBP1-AS1, C21orf62-AS1, ICOSLG, AATBC, MX1, ERG, CYYR1, AP001610.2,
RIPK4, TSPEAR), and 7 genes deviate only in the adjusted run (KCNE1,
OLIG1, RUNX1, AP000282.1, CBR3, ADAMTS1, ATP5PF). 11 genes clear the
magnitude cut in the adjusted model without reaching padj \< 0.01; the
adjusted model has 24 more parameters and 43 fewer samples, so its
adjusted p-values are less powerful as well as adjusted.

| lane_baseline   | lane_adjusted        |   N |
|:----------------|:---------------------|----:|
| DE_high         | DE_high              |   3 |
| DE_high         | Expected_dosage      |   4 |
| DE_high         | Not_DE_outside_noise |   3 |
| DE_low          | DE_low               |   3 |
| DE_low          | Expected_dosage      |   8 |
| DE_low          | Not_DE_outside_noise |   2 |
| Expected_dosage | DE_high              |   3 |
| Expected_dosage | DE_low               |   4 |
| Expected_dosage | Expected_dosage      | 117 |
| Expected_dosage | Not_DE_outside_noise |   6 |
| High_repeats    | High_repeats         |  12 |
| Low_expression  | Low_expression       | 153 |

Lane transitions, baseline to adjusted, all chr21 genes.

| Gene | lane (baseline) | lane (adjusted) | log2FC (baseline) | log2FC (adjusted) | eQTL (baseline) | eQTL (adjusted) |
|:---|:---|:---|---:|---:|:---|:---|
| OLIG2 | DE_low | DE_low | -1.19 | -2.18 | no_cis_eqtl | no_cis_eqtl |
| COL6A1 | DE_low | Not_DE_outside_noise | -0.83 | -0.48 | cis_eqtl | not_evaluated |
| TEKT4P2 | DE_low | Not_DE_outside_noise | -0.64 | -0.53 | no_cis_eqtl | not_evaluated |
| AF165147.1 | DE_low | Expected_dosage | -0.57 | -0.38 | no_cis_eqtl | not_evaluated |
| PDE9A | DE_low | Expected_dosage | -0.55 | -0.15 | cis_eqtl | not_evaluated |
| ZBTB21 | DE_low | Expected_dosage | -0.51 | -0.29 | no_GTEx_data | not_evaluated |
| RBM11 | DE_low | Expected_dosage | -0.47 | 0.21 | no_cis_eqtl | not_evaluated |
| LINC01679 | DE_low | Expected_dosage | -0.45 | -0.27 | cis_eqtl | not_evaluated |
| PAXBP1-AS1 | DE_low | Expected_dosage | -0.45 | -0.09 | no_cis_eqtl | not_evaluated |
| PCBP3 | DE_low | DE_low | -0.43 | -0.50 | cis_eqtl | cis_eqtl |
| BACE2 | DE_low | DE_low | -0.43 | -0.58 | no_cis_eqtl | no_cis_eqtl |
| C21orf62-AS1 | DE_low | Expected_dosage | -0.42 | 0.10 | cis_eqtl | not_evaluated |
| ICOSLG | DE_low | Expected_dosage | -0.42 | -0.14 | cis_eqtl | not_evaluated |
| KCNE1 | Expected_dosage | DE_low | -0.38 | -0.50 | not_evaluated | cis_eqtl |
| OLIG1 | Expected_dosage | DE_low | -0.36 | -0.76 | not_evaluated | cis_eqtl |
| RUNX1 | Expected_dosage | DE_low | -0.31 | -0.43 | not_evaluated | no_GTEx_data |
| AP000282.1 | Expected_dosage | DE_low | -0.27 | -0.70 | not_evaluated | no_GTEx_data |
| CBR3 | Expected_dosage | DE_high | 0.10 | 0.53 | not_evaluated | cis_eqtl |
| ADAMTS1 | Expected_dosage | DE_high | 0.21 | 0.48 | not_evaluated | cis_eqtl |
| ATP5PF | Expected_dosage | DE_high | 0.27 | 0.48 | not_evaluated | no_GTEx_data |
| AATBC | DE_high | Expected_dosage | 0.44 | 0.14 | cis_eqtl | not_evaluated |
| MX1 | DE_high | Not_DE_outside_noise | 0.47 | 0.46 | cis_eqtl | not_evaluated |
| YBEY | DE_high | DE_high | 0.50 | 0.47 | cis_eqtl | cis_eqtl |
| ERG | DE_high | Expected_dosage | 0.50 | -0.29 | no_GTEx_data | not_evaluated |
| COL6A2 | DE_high | DE_high | 0.56 | 0.42 | cis_eqtl | cis_eqtl |
| CYYR1 | DE_high | Expected_dosage | 0.57 | 0.37 | cis_eqtl | not_evaluated |
| AP001610.2 | DE_high | Expected_dosage | 0.64 | 0.39 | no_GTEx_data | not_evaluated |
| RIPK4 | DE_high | Not_DE_outside_noise | 0.64 | 0.45 | cis_eqtl | not_evaluated |
| TSPEAR | DE_high | Not_DE_outside_noise | 0.75 | 0.49 | cis_eqtl | not_evaluated |
| ABCC13 | DE_high | DE_high | 1.45 | 1.01 | cis_eqtl | no_cis_eqtl |

Genes deviating in either run. Ploidy-corrected DESeq2 log2FC; eQTL
columns are empty where a gene was not deviating in that run.

<img src="figures/adjusted_run_comparison_adjusted_vs_baseline.png" alt="One row per gene deviating in either run: ploidy-corrected log2 fold change under the baseline (orange) and adjusted (purple) runs, joined by a line, with a dotted line at zero." width="100%" />

| Gene | largest term | log2 contribution | log2FC (baseline) | log2FC (adjusted) |
|:---|:---|---:|---:|---:|
| OLIG2 | Basophils | 0.448 | -1.19 | -2.18 |
| PDE9A | naïve CD4+ T | -0.415 | -0.55 | -0.15 |
| C21orf62-AS1 | naïve CD4+ T | -0.344 | -0.42 | 0.10 |
| CBR3 | naïve CD4+ T | -0.328 | 0.10 | 0.53 |
| ICOSLG | CD27- B | -0.314 | -0.42 | -0.14 |
| BACE2 | Basophils | 0.277 | -0.43 | -0.58 |
| RBM11 | naïve CD4+ T | -0.274 | -0.47 | 0.21 |
| MX1 | non-classical monocytes | 0.246 | 0.47 | 0.46 |
| COL6A1 | naïve CD4+ T | -0.230 | -0.83 | -0.48 |
| OLIG1 | naïve CD4+ T | 0.224 | -0.36 | -0.76 |
| COL6A2 | CD8+ TEM | 0.217 | 0.56 | 0.42 |
| AP001610.2 | non-classical monocytes | 0.209 | 0.64 | 0.39 |
| TSPEAR | naïve CD4+ T | 0.207 | 0.75 | 0.49 |
| KCNE1 | naïve CD4+ T | 0.193 | -0.38 | -0.50 |
| AF165147.1 | naïve CD4+ T | -0.179 | -0.57 | -0.38 |

Per-covariate attribution of the fold-change shift (S1 panel D, artifact
scale): for each deviating gene, the covariate whose product of
(coefficient) x (T21 minus control mean) is largest. Descriptive;
correlated fractions make single terms less stable than the total.

<img src="figures/adjusted_S1_covariate_evidence.png" alt="Four-panel supplement. A: bar chart of minus log10 p for each candidate covariate against karyotype, coloured by whether it is in the adjusted model. B: number of genes associated with each CyTOF cell fraction within T21 at 5 percent FDR. C: chr21 expression index by karyotype subtype with reference lines at 1 and 1.5; mosaic subjects sit near 1. D: per-covariate contributions to each deviating gene's fold-change shift." width="100%" />
