# Methods spec: Hunter compliance, deviation threshold, and eQTL controls

Decisions from the 2026-08-31 working session. This is the spec that
`docs/superpowers/plans/2026-08-31-threshold-and-eqtl-controls.md`
implements. Companion: `docs/ABSTRACT_DEVIATIONS.md` (what changed after
the stale-input bug was fixed).

---

## Problem 0: the ploidy arm has no library-size normalization

Hunter et al. (2023) list five prerequisites for ploidy-corrected DE. Audit:

| # | Hunter step | Status |
|---|---|---|
| 1 | Minimum read coverage filter (they used 30) | **SATISFIED, as a label.** The absolute count form (total >= 30) excludes nothing: minimum on chr21 protein-coding is 437, median 195,547. The q20 baseMean rule (threshold 25.10) flags 32 of 160 genes into the `Low_expression` lane; Hunter's baseMean < 30 flags 34, of which the 32 are a strict subset. A 2-gene difference. |
| 2 | Mask repeats / remove multi-mapping reads at counting | **SATISFIED at gene level, as a label.** `scripts/archive/process_blacklist.R` overlaps chr21 genes against ENCODE `hg38-blacklist.v2.bed.gz` -> 10 genes, combined with 9 `KNOWN_REPEAT_GENES` into the `high_repeat` flag and the `High_repeats` lane. Read-level masking is unavailable because counts arrive precomputed from Synapse - state as a manuscript limitation. |
| 3 | Remove chr21 for size factor calculation | **RAW arm yes** (`01:175`), **NORM arm no** - see below. |
| 4 | For noisy samples, increase depth or replication | N/A at n = 399; residue is DESeq2 nulling p-values for count outliers (MX1). |
| 5 | Adjust the null / normalize counts by ploidy | **YES** in intent (`01:134-144`), but misapplied - see below. |

**The step 3/5 defect.** `scripts/01_deseq2_analysis.R:210` reads:

```r
normalizationFactors(dds_norm) <- norm_matrix
```

`norm_matrix` contains only 1.0 and 1.5 - line 154 asserts exactly that. In
DESeq2 `normalizationFactors` are the complete per-gene-per-sample divisors and
**replace** size factors; `DESeq()` skips `estimateSizeFactors()` when they are
already set. So the ploidy-corrected arm ran with **no library-size
normalization at all**, on samples spanning 21.1M to 52.3M reads
(2.47x max/min, CV 12.2%). Step 3 is moot there because no size factors exist
to exclude chr21 from.

**Decision.** Replace with what Hunter et al. specify and what the raw arm
already does:

```r
dds_norm <- estimateSizeFactors(dds_norm, normMatrix = norm_matrix,
                                controlGenes = non_chr21_genes_norm)
```

Add a `cooksCutoff = FALSE` sensitivity arm for step 4 and state the step 2
limitation in the manuscript. Steps 1 and 2 need no code change.

**Design principle: label, do not filter.** Hunter et al. *remove* genes failing
steps 1 and 2. This pipeline deliberately *labels* them - `low_expr` and
`high_repeat` feeding the `Low_expression` and `High_repeats` lanes - so genes
below threshold stay in `chr21_lane_assignments.csv` and remain available for
later eQTL-support analysis. Nothing in this spec may convert a label into a
filter. Where a filter is genuinely needed (estimating the robust null, Problem
3 below), it restricts the *estimation set* only; excluded genes stay in the
reported table carrying `q_outlier = NA`.

**This changes every ploidy-corrected number below.** All figures quoted in the
rest of this spec were measured on the unnormalized arm and must be
re-measured after the fix.

**Where we go beyond Hunter et al.:** the MAD/FDR outlier machinery and the eQTL
controls have no counterpart in the paper. Their family-of-four design got an
effect-size filter for free (small n -> large lfcSE -> only large effects
cleared `padj < 0.01`). At n = 399 that protection is gone, so an explicit one
is required. This is an addition to their method, not a departure.

---

## Problem 1: the cohort-noise null is the wrong reference

**Current behaviour.** Scripts 02 and 04 flag a chr21 gene as deviating when
`|norm_log2FC| / SD(non-chr21 protein-coding norm_log2FC) >= 1.0`, i.e.
`|norm_log2FC| >= 0.401`.

**Why it is wrong.** The ploidy-normalization matrix is 1.5 for chr21 genes
in T21 samples and 1.0 everywhere else, so the correction does not act on
diploid genes. Measured:

| set | mean \|raw_log2FC - norm_log2FC\| |
|---|---|
| non-chr21 (n = 15,385) | 0.0048 |
| chr21 (n = 160) | 0.583 (= log2(1.5)) |

A diploid gene's "ploidy-normalized" log2FC is therefore just its ordinary
T21-vs-Control log2FC. Its spread measures real trans effects of trisomy on
diploid genes plus technical noise - a different quantity from "deviation
from the 1.5x dosage expectation". Using it as the yardstick is a category
error, and it is inflated by genuine DS biology.

**Decision.** Estimate the null from the chr21 genes themselves.

## Problem 2: SD is inflated by the signal being detected

`SD/MAD = 1.49` on the chr21 corrected distribution. Removing the 5 most
extreme genes of 160 moves SD by -18.9% but MAD by -1.1%. Those 5 include
OLIG2, a primary compensation candidate: under an SD-based rule OLIG2 helps
set the bar it must then clear.

**Decision.** Use MAD (R's `mad()`, already scaled by 1.4826), not SD.

## Problem 3: the scale must be estimated after the expression filter

log2FC variance depends on expression. Estimating the null over all 160
genes lets near-zero-count genes dominate: the top hits by robust z are
KRTAP12-4 (baseMean 1.4), CLIC6 (5.6), KRTAP10-9 (1.3), CBS (16.7),
KRTAP10-2 (1.3).

**Decision.** Apply the low-expression (`LOW_EXPR_QUANT = 0.20`) and
high-repeat filters FIRST, then estimate median and MAD from the survivors
(n = 119), then test only those genes.

## Problem 4: k = 1 is not a threshold, it is the top third

One SD/MAD-sigma selects ~32% of any normal distribution. Observed at k=1:
36.2% of chr21 genes flagged vs 31.7% expected under the null - a 1.1x
enrichment. Enrichment over the normal null by k: 1.5 -> 1.7x, 2.0 -> 3.1x,
2.5 -> 7.3x, 3.0 -> 19x.

**k can be set statistically.** After trimming the 10% most extreme genes the
chr21 bulk is normal (Shapiro-Wilk p = 0.17; the untrimmed set gives
p = 6e-6). A normal bulk with non-normal tails licenses treating this as a
multiple-testing problem.

**Decision.** Replace the fixed `MAGNITUDE_THRESHOLD` with an FDR criterion:

```
z_i = (norm_log2FC_i - median) / MAD        # median, MAD from eligible genes
p_i = 2 * pnorm(-|z_i|)
q_i = p.adjust(p, "BH")                     # over eligible genes only
```

Primary threshold `OUTLIER_FDR = 0.10`. On current data this yields 5 genes
at an effective k of 2.91 (FDR<0.05 gives 4 genes at k = 3.42).

**Composite qualification rule.** A gene is a deviation call when all three
hold:

1. passes the low-expression and high-repeat filters;
2. `norm_padj < ALPHA` - the deviation is real, not sampling noise;
3. `q_outlier < OUTLIER_FDR` - the deviation is unusually large.

**Reporting.** Always emit a sensitivity table across k (and across
OUTLIER_FDR), never a single count. The lower-than-expected set is
threshold-driven (4 to 18 genes depending on the rule) while the
higher-than-expected set is stable at 6; this asymmetry must be reported.

**Scope limit to state in the paper.** A chr21-internal null can only detect
gene-specific deviation. A uniform chromosome-wide compensation effect would
be absorbed into the median and become invisible.

## Problem 5: MX1's p-value is nulled by an artifact

MX1: baseMean 9,017, norm_log2FC 0.819, robust z 3.73 (second-strongest
signal on the chromosome), but DESeq2 returns `pvalue = NA` and `padj = NA`
while still reporting `stat = 5.63`. It is one of 3 chr21 genes and 217 genes
genome-wide with a nulled p-value - the signature of Cook's-distance outlier
filtering. Under the pre-fix table it had `norm_padj = 1.18e-04`, which is
why the abstract lists it as a positive control.

**Decision.** Diagnose the Cook's flag, then report a `cooksCutoff = FALSE`
sensitivity arm alongside the default. Do not silently reinstate MX1.

## Problem 6: the eQTL "explained" call has no null

`eqtl_lane = explained` requires >= 1 cis variant with (a) GTEx slope sign
matching the deviation sign and (b) within-T21 slope sign matching GTEx at
`t21_p < 0.05`. Three defects:

1. **No multiplicity control.** N cis variants per gene ranges 21 to 1,083,
   with no correction. A gene with 1,083 variants expects ~54 hits at
   p < 0.05 by chance.
2. **Criterion (b) is nearly free.** Measured GTEx-vs-T21 slope agreement is
   93.9% over all variants and 99.7% over T21-significant ones. Almost any
   variant reaching p<0.05 satisfies it, so the call rests entirely on
   criterion (a) - a ~1-in-N coin flip over correlated variants.
3. **N confounds the call.** Median n_cis_total is 107 for "explained" genes
   vs 36 for the one "unexplained" gene (BACE2). RBM11 is called explained on
   a single supporting variant of 83, where chance predicts ~4.

**Decision - controls to add.** All five:

| # | Type | Control | Expected if the criterion is sound |
|---|---|---|---|
| C1 | negative | Genotype permutation: shuffle subject labels, refit | "explained" rate collapses toward 0 |
| C2 | negative | Direction flip: negate `norm_log2FC`, recount | "explained" rate drops substantially |
| C3 | negative | Expected_dosage genes through the identical test | low "explained" rate |
| C4 | positive | GTEx vs within-T21 slope concordance | **DONE, PASSES**: Spearman 0.841, 93.9% sign agreement |
| C5 | positive | Strong chr21 eGenes outside the DE set | high "explained" rate |

**Decision - criterion fixes.**

- **Not** the lead variant. In LD the top-associated variant is frequently a
  tag, not the causal one; substituting it trades one arbitrary rule for
  another.
- **Gene-level permutation p-values** (the GTEx/FastQTL eGene procedure):
  permute sample labels, take min p per permutation, build the null of the
  best variant. Handles LD correlation and variant multiplicity together.
  `p_gene = (1 + #{min_p_perm <= min_p_obs}) / (B + 1)`, B = 1000.
- **BH-FDR across genes** on those gene-level p-values.
- **Quantitative prediction**, not sign agreement: does the eQTL effect
  size predict the *magnitude* of the deviation? Fit
  `observed_deviation ~ predicted_deviation` across genes and report slope,
  R^2, and per-gene residuals. Predicted deviation for gene g uses the
  within-T21 slope and the departure of mean alt dosage from its expectation
  under trisomy.

## Problem 7: there is no Table 1

The manuscript has no cohort characteristics table. One is required, and it
must be defined on the **analysis cohort** - subjects with both whole-blood
RNA-seq and chr21 WGS - not on the expression cohort.

**Decision.** The analysis cohort is **asymmetric**: T21 require RNA-seq **and**
chr21 WGS (302 of 304); Controls require RNA-seq **only** (95 of 95), for a
total of **397**. Genotypes are used solely for the within-T21 dosage
regressions, which controls never enter, so requiring WGS of a control would
discard 89 of 95 (94%) for no analytic gain. The 2 ungenotyped T21 are dropped
so the DE and eQTL analyses cover the same people - previously DE ran on 304 and
eQTL on 302, forcing "302 of 304" phrasing.

Report three columns (T21, Control, overall) plus an explicit row for the T21
excluded for lacking WGS. The roster is built in script 00 from the PASS file
**headers only**, so it is available before any genotype processing and script
01 can subset DE to it.

**Join key.** `subject_id = sub("[A-Z][0-9]*$", "", LabID)`, matching script 02.
Joining on `RecordID` yields zero matches - verified.

- Continuous variables (`Age_at_visit`, `BMI`) as median [IQR] with a Wilcoxon
  rank-sum comparison: both are right-skewed and the control arm is tiny.
- Categorical variables (`Sex`, `Sample_source`, `Event_name`) as n (%) with
  Fisher exact, not chi-square - the genotyped control arm is n = 6, where the
  large-sample approximation fails.
- Optional comorbidity block from `data/P4C_Comorbidity_020921.tsv` (long
  format; one comment line precedes the header row
  `RecordID, Condition, HasCondition, Age.group, min_Age, max_Age`), top 5
  conditions by prevalence.

**Fields that do not exist.** `P4C_metadata_021921_Costello.txt` carries only
`RecordID, Sex, Karyotype, Event_name, LabID, Age_at_visit, BMI,
Sample_source`. There is no race or ethnicity field; record its absence as a
limitation rather than omitting it silently.

## Known data limitation

`data/chr21_ctrl_PASS.csv` contains 14 sample columns, of which **6** overlap
the expression cohort - not 95. There is effectively no internal diploid arm
for eQTL work; GTEx is the only diploid reference. The abstract's "302 people
with T21 ... for comparison" must not imply an internal control group.

Analysis cohort sizes to assert in code:

| arm | RNA-seq | has WGS | **in analysis cohort** |
|---|---|---|---|
| T21 | 304 | 302 | **302** (WGS required) |
| Control | 95 | 6 | **95** (WGS not required) |
| total | 399 | 308 | **397** |

The 6 genotyped controls remain far too few for an internal diploid eQTL arm;
GTEx stays the only diploid reference.

## Out of scope for this plan

- The `normalizationFactors(dds) <- mat` vs
  `estimateSizeFactors(dds, normMatrix = ...)` question in script 01.
- Fine-mapping (SuSiE credible sets) for causal variants.
- Re-rendering the SankeyMATIC diagram (manual browser step).
