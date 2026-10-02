T21-eQTL results summary
================
2026-10-02

Rendered at commit `b10d3d3` (the commit checked out when this document
was last knitted; `results/` is gitignored, so the outputs themselves
are not versioned). This document reads:

- `results/runs/baseline/tables/`: `chr21_lane_assignments.csv` (script
  04\) - every number and gene list in the baseline sections;
  `eqtl_controls_summary.csv` (script 03);
  `positive_control_matching.csv` (script 02);
  `eqtl_best_variant_effects.csv` (script 11).
- `config/runs/baseline.R` - the control thresholds quoted in the text.
- `results/runs/adjusted/`: `processed/cohort_roster.csv` (script 01);
  `tables/chr21_lane_assignments.csv` (script 04);
  `tables/run_comparison_adjusted_vs_baseline.csv` (script 10).
- `docs/figures/` - the PNGs embedded under Figures and in the
  adjusted-run section (copies of the outputs of scripts 06, 07, 10, 11,
  14 and S1), and `Sankey.png`, the SankeyMATIC render of
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
and **3** have no GTEx cis coverage. A cis-eQTL call is a detection
result (q_gene_bh \< 0.05), not an “explained by eQTL” claim.

## Deviating genes

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

In the baseline run, 0 of 100 decoy tests are detected at q \< 0.05 and
10 reach nominal p \< 0.05; 11 of 13 matched positive controls are
detected, beside 14 of 20 tested deviating genes. Both control sets are
drawn in the effect-size figure, and the detection rates in panel D of
the T21-versus-GTEx figure.

## Figures

Copies of the pipeline figures, tracked under `docs/figures/`. To
refresh them, re-run scripts 06, 07, 11, 15, 16, 14 and S1 for each run
(and re-export the lane flow from SankeyMATIC using
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

### Allelic effect in T21 versus GTEx, with the test’s power and base rate (scripts 14, 15, 16):

Panels A and B ask whether each panel’s allele has the same per-copy
effect in T21, fitted at ploidy 3, as in euploid GTEx whole blood, the
assumption behind reading a GTEx eQTL into the trisomic cohort. Panel C
is the power of the gene-level test at each tested deviating gene’s own
GTEx effect size, by spike-in. Panel D is how often the same test
detects a cis-eQTL among the deviating genes, the matched positive
controls and the Expected-dosage genes.

<img src="./figures/baseline_t21_vs_gtex_allelic_effect.png" alt="Four-panel figure. A: scatter of T21 log2 allelic fold change against GTEx whole-blood log2 allelic fold change, per copy of each dosage panel's allele, with 95 percent intervals, for DE high genes (red), DE low genes (blue), positive controls (green triangles) and negative-control decoys (grey diamonds, at GTEx zero); points lie along the dashed identity line. B: the T21-minus-GTEx difference with 95 percent intervals, one row per gene in four groups (DE high, DE low, positive control, negative control), with asterisks where the difference is significant at FDR 0.05. C: spike-in power curves of the gene-level test against the spiked log2 allelic fold change, one curve per tested deviating gene, with a point at the gene's GTEx effect size, filled if detected in T21 and open if not, and a dashed line at 80 percent power. D: fraction of genes with a detected cis-eQTL, with Wilson 95 percent intervals, for deviating genes, matched positive controls, Expected-dosage genes, and Expected-dosage genes in the deviating genes' expression range, with a dashed line at the deviating-gene rate." width="100%" />

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

Of the 23 genes deviating in the baseline, 6 still deviate after
adjustment (OLIG2, PCBP3, BACE2, YBEY, COL6A2, ABCC13), 17 do not
(COL6A1, TEKT4P2, AF165147.1, PDE9A, ZBTB21, RBM11, LINC01679,
PAXBP1-AS1, C21orf62-AS1, ICOSLG, AATBC, MX1, ERG, CYYR1, AP001610.2,
RIPK4, TSPEAR), and 7 genes deviate only in the adjusted run (KCNE1,
OLIG1, RUNX1, AP000282.1, CBR3, ADAMTS1, ATP5PF). 11 genes clear the
magnitude cut in the adjusted model without reaching padj \< 0.01; the
adjusted model has 24 more parameters and 43 fewer samples, so its
adjusted p-values are less powerful as well as adjusted.

<img src="figures/adjusted_run_comparison_adjusted_vs_baseline.png" alt="One row per gene deviating in either run: ploidy-corrected log2 fold change under the baseline (orange) and adjusted (purple) runs, joined by a line, with a dotted line at zero." width="100%" />

<img src="figures/adjusted_S1_covariate_evidence.png" alt="Four-panel supplement. A: bar chart of minus log10 p for each candidate covariate against karyotype, coloured by whether it is in the adjusted model. B: number of genes associated with each CyTOF cell fraction within T21 at 5 percent FDR. C: chr21 expression index by karyotype subtype with reference lines at 1 and 1.5; mosaic subjects sit near 1. D: per-covariate contributions to each deviating gene's fold-change shift." width="100%" />
