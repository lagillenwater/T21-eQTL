# T21 Chromosome 21 Study — Plain-Language Walkthrough

Internal working document. Stays local; not for outside distribution.
Written 2026-09-11 on branch `lncrna-biotype`.

**Primary analysis: the covariate- and composition-adjusted cohort
(`adjusted` run).** Age, weight, and blood-cell-type mix differ between
the Down syndrome and Control groups in this cohort, and blood-cell mix
alone tracks thousands of genes (Supplement A). Reporting the deviation
without accounting for that risks mistaking one of those differences for
a chromosome-21-dosage effect. So the numbers in Steps 1–9 below are all
from the run that statistically adjusts for age, sex, BMI, sample source,
and 19 measured blood-cell-type fractions, and excludes subjects with
mosaic Down syndrome (whose cells don't all carry the extra chromosome).
The unadjusted, Hunter-style analysis on the full cohort is kept and
reported in Supplement B, both for transparency and because it shows what
the adjustment changes.

## The question

People with Down syndrome carry three copies of chromosome 21 instead of
two. Because of that extra copy, genes on chromosome 21 are expected to
make about 1.5 times as much RNA as they do in people with the usual two
copies. Most chromosome 21 genes follow that pattern. Some do not — they
make noticeably more or less RNA than the extra copy alone would predict,
even after accounting for the other ways the two groups differ.

This project asks: for the genes that break the 1.5x pattern, how many
sit near a common DNA variant that is already known to control how much
RNA that gene makes (called a cis-eQTL, cataloged by the public GTEx
project in whole blood)? And how many do not, leaving them as open
candidates for a different kind of explanation — one involving how the
extra chromosome itself changes gene control, rather than a DNA variant
that was already there. Finding a nearby variant only describes a gene;
it does not prove the variant causes the deviation.

---

## Step 1 — Building the study group

**Scripts:** `00_preprocess_data.R`, run configuration in `config/runs/adjusted.R`

**Method.** We combined blood RNA sequencing data with a spreadsheet of
each person's age, sex, weight, and other clinical details, matched people
between the two, and recorded each person's karyotype subtype (typical
trisomy, mosaic, translocation, or unspecified) from a separate genetic
reference. That gave a base group of 397 people with both RNA and
chromosome 21 DNA data (302 Down syndrome, 95 Control). For our primary
analysis, we then required every person to also have a complete set of
the factors we adjust for — age, sex, BMI, sample source, and a matching
blood-cell-type measurement — and excluded anyone with mosaic Down
syndrome, since their cells don't uniformly carry the extra chromosome and
the 1.5x expectation doesn't cleanly apply to them.

**Result.** The primary study group is **354 people: 274 with Down
syndrome, 80 without.** Compared with the base group of 397, this dropped
22 people missing a blood-cell-type measurement, 12 missing a BMI
measurement, 9 with mosaic Down syndrome, and 2 without chromosome 21 DNA
data. In this group, BMI, blood sample source, and study visit still
differ significantly by karyotype (each p < 0.001–0.005), which is exactly
why we adjust for them rather than ignore them; age and sex do not differ
significantly. Several health conditions remain far more common in the
Down syndrome group (underactive thyroid, sleep apnea, obesity, autoimmune
skin conditions, congenital heart defects — all p < 0.0001), as expected.

**Table:** `results/runs/adjusted/tables/table1_analysis_cohort.md` — full
characteristics table for the primary 354-person cohort. The full,
pre-adjustment 397-person cohort is described in Supplement A. No figure
for this step.

---

## Step 2 — Correcting for the extra chromosome and for other differences between groups

**Script:** `01_deseq2_analysis.R`

**Method.** Comparing gene activity between two groups normally assumes
"no difference" means a 1-to-1 ratio. That's wrong for chromosome 21 in
Down syndrome, where a perfectly normal gene will still read as "elevated"
simply because there's 50% more DNA to read from. We corrected this by
telling the statistical model to expect a 1.5x ratio on chromosome 21
specifically (in Down syndrome samples only) and a 1-to-1 ratio elsewhere,
and by excluding chromosome 21 from the calculation that makes samples
comparable to each other in the first place. On top of that
ploidy correction, the primary model's statistical design also includes
age, sex, BMI, sample source, and the 19 blood-cell-type fractions, so the
karyotype effect it reports is estimated net of those factors. The
genetic-variant test in Steps 3–4 uses a matching version of gene activity
with those same covariate effects mathematically removed.

**Result.** The correction works as intended: before it, nearly every
chromosome 21 gene sits above the 1.5x line and reads as statistically
significant; after it, chromosome 21 genes are centered on "no difference
from 1.5x" and mixed in with the rest of the genome (quantified in Step
7). A PCA plot of all samples (a standard check of how similar samples are
to each other overall) shows Down syndrome and Control samples
overlapping rather than forming two separate clusters, confirming that
karyotype alone isn't the dominant source of variation in blood gene
activity — supporting the decision to adjust for the other factors that
are.

**Figures:**
- `results/runs/adjusted/figures/ma_plot.pdf` — four panels: raw vs.
  corrected, all genes vs. chromosome 21 highlighted. Chromosome 21 genes
  (green) sit above the 1.5x line before correction and collapse to match
  the rest of the genome after.
- `results/runs/adjusted/figures/pca_plot.pdf` — samples colored by
  karyotype; no clean separation.
- `results/runs/adjusted/figures/dispersion_plot.pdf` — standard
  model-fit quality check.

---

## Step 3 — Picking out the genes that break the pattern, and pulling their genetic data

**Script:** `02_filter_genotypes.R`

**Method.** From all chromosome 21 genes, we kept the three gene types
GTEx tests for cis-eQTLs (protein-coding, lncRNA, pseudogene — 318 genes
total) and set aside genes with too little RNA to measure reliably (fewer
than 30 average reads) or genes in repetitive DNA that's hard to measure
accurately. Among the rest, in the adjusted model, a gene counts as
**deviating** if its difference from the 1.5x line is statistically
convincing (adjusted p-value below 0.01) and at least 33% away from the
line. Genes 50% or more away get the strongest label ("tier 1"); genes
between 33% and 50% away get a preliminary label ("tier 2"). For every
deviating gene, we pulled every DNA variant within about 1 million DNA
letters of it that GTEx had already flagged as at least loosely linked to
its activity in whole blood, and looked up each genotyped Down syndrome
subject's DNA at those spots.

**Result.** **13 of the 318 genes deviate** in the primary analysis: 6
make more RNA than expected and 7 make less. 4 are tier 1 (strong), 9 are
tier 2 (near-threshold). This is fewer than the 23 genes that deviate in
the unadjusted analysis (Supplement B) — accounting for age, weight, and
blood-cell mix removes some genes from the deviating set and adds others,
discussed in Step 5 and Supplement B.

**Table:** `data/processed/eqtl_supported_genes.csv`,
`eqtl_target_variants.csv`. No figure for this step; the deviating genes
are visualized in Steps 6 and 8.

---

## Step 4 — Testing whether a nearby genetic variant tracks each gene's RNA level

**Script:** `03_t21_dosage_boxplots.R`

**Method.** For each candidate variant near a deviating gene, we checked
whether Down syndrome subjects with more copies of that variant (0 to 3)
tended to have higher or lower RNA levels for that gene, on the
covariate-adjusted expression measure from Step 2. For each gene we kept
its single strongest result, then judged how strong that result really
was by shuffling the RNA measurements across people 1,000 times and
rerunning the same search each time — this shows how often a result this
strong would appear by chance. A gene counts as having a **detected
cis-eQTL** if fewer than 5% of the 1,000 shuffles beat its real result,
after adjusting for testing multiple genes at once (the same "eGene"
method GTEx itself uses). We also ran two trustworthiness checks through
the identical procedure: a **negative control** (each deviating gene
tested against variants near a different gene at least 5 million DNA
letters away, which should rarely look like a hit) and a **positive
control** (the 10 chromosome 21 genes with the strongest known GTEx
signal, which should mostly be caught).

**Result.** Of the 13 deviating genes, 3 have no GTEx-tested variant
nearby and can't be tested (ATP5PF, RUNX1, AP000282.1). Of the 10
testable genes, **7 have a detected cis-eQTL and 3 do not** (ABCC13,
BACE2, OLIG2). The negative control found nothing (0 of 10), at or below
the 5% false-detection rate expected from chance. The positive control
caught all 10 of 10 known genes. Typical effect sizes on the deviating
genes' own variants were modest (about 0.13–0.20 log2-CPM per extra copy
of the variant) and only somewhat larger than the negative-control decoys
(about 0.12) but clearly smaller than the positive controls (about 0.48).

A detected cis-eQTL describes the gene — its RNA level in Down syndrome
blood tracks a common DNA variant nearby. It does not, by itself, prove
that variant is *why* the gene deviates from 1.5x.

**Tables:** `results/runs/adjusted/tables/eqtl_gene_level_perm.csv`,
`eqtl_control_negative.csv`, `eqtl_control_positive.csv`,
`eqtl_controls_summary.csv`. Figures for this step are in Step 9.

---

## Step 5 — Sorting every chromosome 21 gene into a final category

**Script:** `04_chr21_lane_assignment.R` — **this is the main result
table.**

**Result.**

| Category | Genes | What it means |
|---|---|---|
| Expected dosage | 129 | Matches the 1.5x expectation |
| Low expression | 153 | Too little RNA measured to assess |
| High repeats | 12 | In repetitive DNA, not reliably measurable |
| Near-threshold, not significant | 11 | Off the 1.5x line by 33%+ but short of the p < 0.01 bar — see note below |
| Deviates, higher (DE_high) | 6 | See below |
| Deviates, lower (DE_low) | 7 | See below |

Of the 13 deviating genes:

| Gene | Direction | Strength | Genetic signal (cis-eQTL) |
|---|---|---|---|
| ABCC13 | Higher | Strong (tier 1) | Tested, not detected |
| CBR3 | Higher | Near-threshold (tier 2) | Detected |
| ADAMTS1 | Higher | Near-threshold (tier 2) | Detected |
| YBEY | Higher | Near-threshold (tier 2) | Detected |
| COL6A2 | Higher | Near-threshold (tier 2) | Detected |
| ATP5PF | Higher | Near-threshold (tier 2) | Not testable — no GTEx data |
| OLIG2 | Lower | Strong (tier 1) | Tested, not detected |
| OLIG1 | Lower | Strong (tier 1) | Detected |
| AP000282.1 | Lower | Strong (tier 1) | Not testable — no GTEx data |
| BACE2 | Lower | Near-threshold (tier 2) | Tested, not detected |
| KCNE1 | Lower | Near-threshold (tier 2) | Detected |
| PCBP3 | Lower | Near-threshold (tier 2) | Detected |
| RUNX1 | Lower | Near-threshold (tier 2) | Not testable — no GTEx data |

**Six genes remain open candidates** for a mechanism other than a common
cis-eQTL: 3 were tested and showed no genetic signal (ABCC13, BACE2,
OLIG2), and 3 could not be tested for lack of GTEx data (ATP5PF, RUNX1,
AP000282.1). This does not rule out a cis-eQTL for these genes; it means
the pipeline has not found one.

**On the 11 "near-threshold, not significant" genes:** these clear the
33% magnitude bar but not the p < 0.01 bar in this smaller, more heavily
adjusted model (354 people and 24 more statistical terms, versus 397
people and no extra terms for the unadjusted model). Several of them do
deviate under the unadjusted analysis (Supplement B); their disappearance
here is at least partly a loss of statistical power, not necessarily
evidence the underlying effect is gone.

**Table:** `results/runs/adjusted/tables/chr21_lane_assignments.csv` (the
canonical per-gene table), `chr21_lane_summary.csv`. No figure for this
step; see Steps 6, 8, and 9.

---

## Step 6 — Drawing the flow diagram

**Script:** `05_sankeymatic_export.R`, rendered at sankeymatic.com

**Method.** We turned the Step 5 category counts into a flow (Sankey)
diagram: all 318 genes on the left, splitting into "Expected dosage," "Not
assessable," "Near-threshold, not significant," and "Outside expectation"
(which splits into higher/lower, then into the three cis-eQTL outcomes).

**Result.** One figure showing the whole classification at a glance.

**Figure:** the input file
`results/runs/adjusted/tables/chr21_lane_sankeymatic_input.txt` is ready,
but **has not yet been pasted into sankeymatic.com and rendered as an
image** — that's a manual step. The rendered baseline version
(`docs/figures/Sankey.png`) exists for the unadjusted analysis; an
adjusted-cohort version still needs to be generated the same way before
this can go in a manuscript.

---

## Step 7 — Checking that the correction actually worked

**Script:** `06_chr21_distribution_panel.R`

**Method.** We plotted the distribution of gene-activity differences
(Down syndrome vs. Control) for all 318 targeted chromosome 21 genes,
before and after the Step 2 correction, and compared it against
chromosome 22 — a chromosome with the usual two copies in everyone, which
the correction should leave untouched.

**Result.** Before correction, the typical chromosome 21 gene sits almost
exactly at the 1.5x line. After correction, that typical value shifts down
to almost exactly zero — the shift (0.585 on the log2 scale) matches the
1.5x expectation almost perfectly. Chromosome 22's distribution does not
move (shift effectively zero), confirming the correction is specific to
chromosome 21.

**Figure:** `results/runs/adjusted/figures/ploidy_correction_distributions.pdf`
**Table:** `results/runs/adjusted/tables/ploidy_correction_distribution_stats.csv`

---

## Step 8 — Volcano plots

**Script:** `07_three_panel_figure.R`

**Method.** Volcano plots show, for every gene, how big its expression
difference is (x-axis) against how statistically confident that
difference is (y-axis). One plot compares all target-type genes
genome-wide before/after correction with chromosome 21 highlighted; the
other shows only chromosome 21, with the 13 deviating genes labeled.

**Result.** In the first plot, the uncorrected chromosome 21 genes form a
tight, almost entirely "significant" cluster near the 1.5x line; after
correction that cluster centers on zero and blends into the rest of the
genome. The second plot shows exactly which 13 genes deviate and in which
direction (red = higher, blue = lower, bold = tier 1).

**Figures:**
- `results/runs/adjusted/figures/volcano_all_genes.pdf`
- `results/runs/adjusted/figures/volcano_chr21.pdf`

---

## Step 9 — Figures for the genetic-signal test

**Script:** `11_eqtl_figures.R`

**Method.** Four figures from the Step 4 results: (1) for each deviating
gene, RNA level against the number of copies of its strongest-linked
variant, genes higher than expected in panel A and lower than expected in
panel B; (1b) the same layout as a supplement, with the positive controls
in panel A and the negative controls in panel B; (2) a scatter of effect
size against statistical confidence for every gene and its controls; (3) a
map of chromosome 21 marking each deviating gene and whether it had a
detected genetic signal.

**Result.** The dosage panels make the Step 4 numbers visible gene by
gene. The effect-size scatter shows the positive controls spread widely
with large effects, the deviating genes with more modest effects, and the
negative controls sitting lower — though with less separation from the
deviating genes than in the unadjusted analysis (Supplement B), consistent
with the smaller, more heavily adjusted sample carrying less statistical
power. The chromosome map shows the 13 deviating genes spread across the
chromosome rather than clustered in one region.

**Figures:**
- `results/runs/adjusted/figures/eqtl_dosage_panels.pdf` (A: DE high, B: DE low)
- `results/runs/adjusted/figures/eqtl_dosage_controls.pdf` (supplement; A: positive controls, B: negative controls)
- `results/runs/adjusted/figures/eqtl_effect_sizes.pdf`
- `results/runs/adjusted/figures/chr21_deviating_map.pdf`

**Table:** `results/runs/adjusted/tables/eqtl_best_variant_effects.csv`

---

## Supplement A — Why we adjust: evidence the groups differ beyond karyotype

**Script:** `S1_covariate_evidence.R`

**Method.** Before trusting that a gene's deviation reflects chromosome 21
dosage rather than age, weight, or blood-cell mix, we checked, on the full
397-person cohort: (A) do age, weight, and other clinical factors really
differ by karyotype; (B) does directly measured blood-cell-type mix
associate with gene activity broadly; (C) do the 9 mosaic Down syndrome
subjects — who don't carry the extra chromosome in every cell — show a
smaller chromosome 21 effect, as a sanity check that the dosage
measurement itself is real; (D) for each gene affected by adjustment,
which single covariate contributed most to the shift.

**Result.** Yes to all three checks: several clinical factors differ by
karyotype (full-cohort Table 1, below), cell-type mix associates with
thousands of genes, and mosaic subjects show a chromosome 21 activity
level of about 1.07x rather than 1.35x for full trisomy — the smaller
effect expected if only some of their cells carry the extra chromosome.
This is the direct evidence behind the decision to adjust in Steps 1–9.
Panel D traces, gene by gene, which covariate (often a specific blood cell
type, e.g. naïve CD4+ T cells or basophils) contributed most to that
gene's shift between the unadjusted and adjusted results.

**Figure:** `results/runs/adjusted/figures/S1_covariate_evidence.pdf`
(all four panels; panel D only exists in the adjusted run's copy since it
compares against the unadjusted result)

**Table (full, pre-adjustment cohort):** `results/tables/table1_analysis_cohort.md`
— 397 people (302 Down syndrome, 95 Control), before any covariate or
mosaic exclusion. Age (p = 0.035) and BMI (p < 0.0001) differ by karyotype
here too; sex does not.

---

## Supplement B — Comparison to the unadjusted (Hunter-style) analysis

**Scripts:** the full pipeline under `config/runs/baseline.R`, then `10_compare_runs.R`

**Method.** We also ran the identical pipeline with no covariate
adjustment and no mosaic exclusion — the original Hunter et al. approach —
on the full 397-person cohort, and compared its results gene by gene
against the adjusted, primary analysis.

**Result.** The unadjusted analysis finds **23 deviating genes** (10
higher, 13 lower) instead of 13, and detects a cis-eQTL for 14 of the 20
testable ones instead of 7 of 10. Of the 13 genes that deviate in the
primary (adjusted) analysis, 6 also deviate unadjusted (OLIG2, PCBP3,
BACE2, YBEY, COL6A2, ABCC13); the other 7 only deviate once age, weight,
and cell mix are accounted for (KCNE1, OLIG1, RUNX1, AP000282.1, CBR3,
ADAMTS1, ATP5PF). Of the 23 unadjusted-only or shared deviating genes, 17
drop out once adjusted. As Step 5 notes, some of that movement reflects
reduced statistical power in the smaller, more heavily adjusted model
rather than a real disappearance of the effect — 11 genes still show a
similar-sized shift in the adjusted model without reaching its stricter
significance bar. This comparison is what justifies reporting the
adjusted analysis as primary while keeping the unadjusted one visible: it
shows concretely which genes' status depends on accounting for these
other factors.

**Figures:**
- `results/runs/adjusted/figures/run_comparison_adjusted_vs_baseline.pdf`
- Full unadjusted-analysis figure set (all of Steps 2–9, run without
  adjustment): `results/runs/baseline/figures/`

**Tables:** `results/runs/adjusted/tables/run_comparison_adjusted_vs_baseline.csv`,
`lane_transitions_adjusted_vs_baseline.csv`,
`results/runs/baseline/tables/chr21_lane_assignments.csv` (full unadjusted
per-gene table)

---

## Supplement C — Confirming the pipeline reproduces itself

**Script:** `audit_baseline_drift.R`

**Method.** We keep a saved copy of a complete earlier pipeline run
(September 4, 2026, unadjusted configuration). Every time the pipeline is
regenerated, an automatic check compares every output table, line by
line, against that saved copy. This is a check on the pipeline's
mechanics, run against the unadjusted configuration because that's the
one with an archived reference copy to check against — it is not itself
about which run is reported as primary.

**Result.** The regenerated results matched the saved copy exactly, with
one intentional exception: one gene, CYYR1, changed from "no genetic
signal found" to "genetic signal found" after we switched how gene
activity is measured for this test (from raw sequencing counts to a
measure adjusted for each person's total sequencing depth) — a deliberate
methods change, not a bug. Every other number in the pipeline is
identical, giving confidence the pipeline is reproducible and that
results only change when something is intentionally changed.

**No figure.** Full comparison log: `results/archive/2026-09-04_flat/drift_audit.md`

---

## Where everything lives

- Main (adjusted, primary) result table:
  `results/runs/adjusted/tables/chr21_lane_assignments.csv`
- Main figures (PDF, publication quality): `results/runs/adjusted/figures/`
- Unadjusted comparison run: `results/runs/baseline/tables/` and
  `results/runs/baseline/figures/`
- Full methods detail and thresholds: `CLAUDE.md` (repository root) — note
  that file currently calls the unadjusted `baseline` run "the primary
  result," which is the pipeline's internal/historical framing, not the
  reporting choice used in this document
- Design decisions and history: `docs/decisions.md`
- Machine-generated numeric summary (denser, technical, both runs):
  `docs/summary.md`

**Open item:** the adjusted-cohort Sankey diagram (Step 6) still needs to
be rendered from `results/runs/adjusted/tables/chr21_lane_sankeymatic_input.txt`.
