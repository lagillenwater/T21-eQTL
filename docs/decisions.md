# Decisions, legacy notes, and gotchas

The README documents the current pipeline only. This file records how it
got there: design decisions and retired approaches, legacy data sources
and terminology, the archived script catalog, and practical gotchas.

## Decision log

### Classification threshold: Hunter's rule; cohort-SD filter retired

The classification rule matches Hunter et al.'s own criterion directly:
`norm_padj < ALPHA_DE` (0.01) AND `abs(norm_log2FC) >= DEVIATION_LFC`
(`log2(1.5)`), applied on the ploidy-corrected scale; that is tier 1. A
labelled tier 2 (`DEVIATION_LFC_T2 = log2(4/3)`, adopted 2026-09-01)
reports near-threshold genes rather than hiding them, and it is the cut
the `sig_lane` split and the script 02 target selection are actually
made at; `DEVIATION_LFC` sets only the `tier` column. An earlier pipeline
instead gated deviations on a cohort-derived noise threshold (1 SD of the
non-chr21 cohort noise; the constant was then named `MAGNITUDE_THRESHOLD`).
That filter was retired: ploidy normalization only acts on chr21 genes, so
a non-chr21 cohort-noise SD is not a valid reference for the chr21 null,
and retiring the filter removed that mismatch. `DEVIATION_LFC_T2` (with
`DEVIATION_LFC` labelling the tiers) is the one dial in this analysis.

### Target biotypes: protein-coding, lncRNA, pseudogene; low-expression floor: Hunter's absolute 30

The chr21 target set is protein-coding, lncRNA and pseudogene genes
(`TARGET_BIOTYPES` in `scripts/lib/biotypes.R`, sourced by scripts 02, 04,
06, 07), because GTEx whole blood tests all three for cis-eQTLs: of the 223
chr21 genes in the GTEx v10 allpairs file, 137 are protein-coding, 61 lncRNA
and 13 pseudogenes in the HTP annotation. An earlier run
restricted to protein-coding (`RESTRICT_TO_PROTEIN_CODING`, now removed),
which left the eQTL-testable lncRNAs and pseudogenes out of the question. Widening the set
is method application, not a new method: the same rules run on more genes.

Widening exposed a dependency in the low-expression flag. It had been the
20th percentile of baseMean within the target set (`LOW_EXPR_QUANT`), which
sat at 25.1 for protein-coding genes but would drop to 4.1 with lncRNA in
the set (chr21 lncRNA median baseMean is 8 against 482 for protein-coding),
admitting genes with a handful of counts. It is now Hunter et al.'s own
absolute coverage floor, `LOW_EXPR_BASEMEAN = 30`, applied identically in
scripts 02 and 04. On the protein-coding set the absolute floor flags 34
genes against the quantile's 32, a strict superset: TFF3 (baseMean 27.5)
and C21orf59-TCP10L (25.2) move from Expected dosage to Low expression, and
no DE call changes.

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
the classification rule: checked on 2026-08-31 (then
`eqtl_negative_controls.csv`, since replaced by the standalone controls
below), it called 15 of 20 genes "explained" as observed, 18 of 20 with
the deviation direction flipped, and 20% with genotypes shuffled, so it
does not discriminate real signal from chance.
The classification rule is the gene-level permutation test
(`scripts/lib/eqtl_fit.R`, run in script 03): for each deviating gene, the
best-variant test statistic is compared against its null distribution
under permutation of genotype-to-expression assignment, giving
`p_gene_perm`; BH-adjusted across deviating genes to `q_gene_bh`, with
`eqtl_lane = cis_eqtl` when `q_gene_bh < FDR_GENE` (0.05).

### Standalone controls for the eQTL test (2026-09-10)

The permutation null inside the gene-level test shuffles expression across
subjects, which is a negative-control operation, but it is the machinery
that produces the p-value rather than a check reported beside it. The
user asked for distinct, standalone positive and negative controls. Both
go through the same runner as the observed set
(`scripts/lib/eqtl_controls.R::gene_level_tests`, script 03), so the
three sets are tested identically:

- **Negative, unlinked variant sets.** Each deviating gene's expression is
  paired with the cis variant set of another tested gene with a TSS at
  least `decoy_min_distance` (5 Mb) away, beyond any LD with its own
  locus; among those the gene with the closest variant count is chosen so
  the decoy carries the same multiplicity (a first version took the
  farthest gene, which handed 16 of 20 genes the same one-variant decoy).
  Real genotypes, the same subjects, the same test; detections should sit
  near the FDR level. `eqtl_control_negative.csv`.
- **Positive, strong GTEx eGenes.** Script 02 selects the
  `n_positive_controls` (10) chr21 genes with the smallest GTEx
  whole-blood nominal p among expressed, non-repeat, non-deviating genes
  and pulls their cis variants alongside the deviating genes
  (`gene_set == "positive_control"`); script 03 tests them on their own
  variants. Most should be detected, or the test lacks power at this
  cohort size. `eqtl_control_positive.csv`.

The old direction-flip and genotype-shuffle checks of the retired
any-variant rule (`eqtl_negative_controls.csv`) were removed; their
numbers are recorded above. Positive-control genes are excluded from
`t21_dosage_per_variant.csv`, `eqtl_gene_level_perm.csv` and the
representative-variant table, so script 04 and the lane table are
unaffected. The observed test itself did not change: same variants, seeds
and BH, now called through the shared runner.

### Figures split by question (2026-09-10)

The 2x2 `Chr21_DEG` volcano coloured deviating genes by direction and eQTL
outcome at once, and the per-gene dosage boxplots lived in script 03. Both
were replaced: script 07 now writes `volcano_all_genes` (chr21 highlighted,
no labels: the ploidy-correction shift) and `volcano_chr21` (deviating genes
labelled, coloured by direction only), and the new script 11 draws the
eQTL stage from existing tables: `eqtl_dosage_panels` (best variant of the
permutation test per gene, with the positive and negative controls in their
own groups), `eqtl_effect_sizes` (a scatter of the within-T21 slope at each
gene's best variant against -log10 of its gene-level permutation q, for
the deviating genes, their decoy sets and the positive controls; earlier
versions plotted q alone, then absolute slopes as a faceted and then a
single forest plot, and were revised at the user's request) and `chr21_deviating_map` (bands at each deviating
gene's TSS coloured by eQTL outcome, labels by direction). Positions come
from the roster TSS, with `data/chr21_gene_positions.csv` (Ensembl GRCh38)
covering deviating genes GTEx does not carry. Nothing statistical moved.

### Dosage panels split from their controls (2026-09-11)

`eqtl_dosage_panels` originally stacked all four groups (DE high, DE low,
positive controls, negative controls) into one figure. At the user's
request it now writes two: `eqtl_dosage_panels` (DE high as panel A, DE
low as panel B) for the main text, and a new `eqtl_dosage_controls`
supplement (positive controls as panel A, negative controls as panel B),
both lettered with `patchwork::plot_annotation(tag_levels = "A")`. Same
underlying data (`panel_variants()` in `scripts/lib/eqtl_figures.R` is
unchanged) and the same per-gene panels; only which groups share a figure
changed.

### eQTL direction stated per minor allele (2026-09-22)

Manuscript issue #3 (greenelab/T21-cis-eqtl-Chr21-regulation-manuscript):
"Align common eQTL minor or major allele with direction in the typical
population." The pipeline codes genotype as ALT-allele dosage (script 02
`alt_dosage`) and GTEx's `slope` is ALT-referenced too, so the two were
internally consistent - but ALT is whichever base differs from the
reference assembly, not the rarer allele. In the baseline cis-variant
universe 19.2% of the retained (variant, gene) pairs have ALT as the
MAJOR allele, and 13 of the 47 variants the dosage panels are drawn on.
A slope "per alt allele" therefore pointed one way for some variants and
the other way for others and could not be compared with a published eQTL
direction.

Every reported direction is now stated per copy of the MINOR allele
(`scripts/lib/alleles.R`). The genotype coding is unchanged; the
alignment is a reflection of the regressor, which flips a slope's sign
and leaves |t|, p, SE and R-squared alone. So `p_gene_perm`, `q_gene_bh`,
`cis_eqtl_detected` and every lane are invariant, and the re-run
reproduced the lane table exactly: 8 DE_high cis_eqtl / 2 no_GTEx_data,
6 DE_low cis_eqtl / 6 no_cis_eqtl / 1 no_GTEx_data, 130 Expected dosage,
12 High repeats, 153 Low expression. Script 11's refit-p check still
matches the stored best-variant p for all 50 panels, on the reflected
dosage.

Which allele is minor is read from two reference populations, both
carried per variant:

- **GTEx whole blood `af`** (the allpairs column), primary: the
  population the cis-eQTL itself was called in.
- **gnomAD v4.1 genomes** (`AF` global and `AF_nfe`), the independent
  check (`scripts/lib/gnomad.R`, `scripts/fetch_gnomad_af.R`). The 7.8 GB
  sites VCF is never downloaded; `Rsamtools::scanTabix` pulls only the
  blocks covering the cis positions over https and the result is cached
  under `data/gnomad/`.
- The HTP T21 cohort's own ALT frequency (alt copies out of three) rides
  along as a third column, so a cohort that does not match either
  reference would be visible.

Over the 5,700 tested baseline variants: gnomAD agrees with GTEx on the
minor allele for 97.3% and the HTP cohort for 98.7%. Every one of the 156
gnomAD disagreements sits at GTEx MAF 0.425 to 0.499 (126 of them above
0.45) - variants where the two alleles are near-equally common and
"minor" is close to a coin flip, so the sign of the reported direction
depends on the population chosen. `eqtl_allele_alignment.csv` carries the
per-variant flags (`minor_concordant`, `htp_concordant`, `maf_tie`) so
those cases are visible rather than buried. One deviating gene's best
variant is affected: BACE2 (GTEx MAF 0.484, ALT minor; gnomAD calls ALT
major).

**Dosage panels are oriented to the deviation-matching allele
(2026-09-22).** Once every direction was stated per minor allele, the
dosage panels became hard to read: inside the DE high block four panels
trended up and four down, because the minor allele's effect direction is
a property of the variant and has no reason to follow the gene's
deviation. Script 11's panels are now drawn on the DEVIATION-MATCHING
allele - the allele whose GTEx effect runs the same way the gene deviates
- so a DE high panel trends up and a DE low panel trends down wherever
the within-T21 fit reproduces GTEx, and a panel that fails to reproduce
it is the one running the wrong way. The orientation is taken from GTEx,
never from the within-T21 slope being plotted; orienting on the data
would force the right trend by construction and show nothing.

The minor allele is plotted exactly when it is itself the deviation-
matching one, so `plot_allele_role` doubles as the with/against-the-
deviation badge (`orient_panels()` in `scripts/lib/eqtl_figures.R`). In
baseline that is 8 of 20 panels, leaving 12 drawn on major-allele dosage
(adjusted: 3 of 10). Every strip still names the minor allele, its MAF
and the canonical GTEx direction, and carries a trend line coloured by
whether the within-T21 fit reproduces that direction - 20 of 20 panels do,
in both runs. Control panels have no deviation to match and stay on
minor-allele dosage; a decoy variant is an eQTL of the decoy gene rather
than of the gene plotted, so it gets no canonical direction. The
effect-size scatter keeps the minor allele as its common cross-gene
reference: orienting it per gene would make the sign of x a restatement
of the block and of the GTEx agreement.

Layout: one row per category, so a block's panels sit side by side and the
figure grows sideways rather than wrapping (at the user's request,
2026-09-22). Width is 3.1 in per panel, the space the four-line strip needs
without clipping, so the adjusted run's 5 + 5 deviating-gene panels give a
16 in figure and baseline's 8 + 12 give 38 in; the controls supplement, with
20 decoy panels in its widest block, reaches 63 in.

The with/against badge is descriptive. T21 subjects carry no excess of
these alleles, so it records a coincidence of direction, not an account
of the deviation - the same framing as `eqtl_lane`.

ALT-referenced columns are kept beside the aligned ones rather than
replaced, so nothing downstream changes meaning silently:
`gtex_slope`/`t21_slope` and `dir_match` stay ALT-referenced, and
`gtex_slope_minor`/`t21_slope_minor` and `dir_match_minor` are the
minor-allele statements. `supportive` compares two ALT-referenced slopes
and is unaffected either way - flipping both signs preserves their match.

### `eqtl_lane` is a detection result, not an "explained by eQTL" claim

`cis_eqtl` means a cis-eQTL is detectable for the gene at
`q_gene_bh < 0.05` in the within-T21 data - not that the eQTL
quantitatively accounts for the observed deviation.

### Descriptive framing: the annotations describe deviations, they do not explain them

Every chr21 gene is placed by its deviation from the trisomy expectation
(ploidy-corrected log2FC against 0). Deviating genes then carry the
cis-eQTL annotation, which is observed rather than inferred and is not a
cause.

The cis-eQTL call has that status for a structural reason. A cis-eQTL
describes between-individual variance at a locus; a lane assignment is a
difference in group means, and common-variant genotype frequencies do not
differ between T21 and control subjects drawn from the same population. A
cis-eQTL therefore cannot produce a mean deviation from 1.5x. "cis-eQTL
detected" says the gene's expression in T21 blood tracks common cis
variation (a property of the gene), not that the variant accounts for the
lane. Hunter et al. could invoke eQTLs because their single T21 individual
had a single genotype; with 302 genotyped subjects that argument does not
carry over.

### Co-expression neighborhood check retired (2026-09-10)

From 2026-08 to 2026-09-10 every deviating gene also carried a
co-expression neighborhood check: its 20 most correlated non-chr21 genes
were chosen in controls, their median T21-vs-Control log2FC was compared
with a correlation-matched null of 300 random modules, and the gene was
labelled by whether its neighbors shifted with it (`neighborhood`,
`partner_z`, `gene_minus_partner_lfc`; script 04 Step 3b,
`scripts/lib/neighborhood.R`, `chr21_neighborhood_shift.csv`). Scripts 08
and 09 then reported the within-T21 R-squared of each deviating gene on its
neighborhood score, beside its best cis variant and against background
genes and controls. The labels were renamed twice (PROGRAM / MIXED /
GENE-SPECIFIC, then SHARED / PARTLY-SHARED / NOT-SHARED, then
neighbors_shift / neighbors_shift_less / neighbors_no_shift) because each
earlier name claimed a cause the test does not observe.

The check, and scripts 08 and 09 with it, were retired on 2026-09-10. In
whole blood a gene's co-expression partners are largely a cell-type
signature, so the check was a proxy for composition, and the adjusted run
now uses measured CyTOF cell fractions directly as covariates
(`config/runs/adjusted.R`, `S1_covariate_evidence.R`). What the proxy
could add was hard to read: a high neighborhood R-squared is the default
for any expressed whole-blood gene (median 0.63 across 1000 random
non-chr21 genes), a low one could mean decoupling in T21 or a poorly chosen
partner set, and on the adjusted artifact, with the cell fractions already
removed, the composition reading no longer applied. Lane classification,
`eqtl_lane` and every DESeq2 table are unchanged by the removal; the
archived 2026-09-04 outputs keep the neighborhood tables for the record.
The matched-null design is described below for that record.

### Named runs, covariate adjustment as one artifact, mosaic exclusion (2026-09-10)

The pipeline had no configuration surface: every path was a literal, the
DESeq2 design was fixed to `~ karyotype`, and thresholds were repeated in
three scripts. Three findings on 2026-09-10 required a covariate-adjusted
run: age, BMI and sample source differ by karyotype (BMI median 27.2 vs
24.3; NDSC2018 source 92 T21 vs 7 control; Event_name also imbalanced),
CyTOF cell composition differs by karyotype and associates with thousands
of genes, and 9 of 302 T21 subjects are mosaic (INCLUDE MONDO codes), with
a chr21 expression index near 1.07 against 1.35 for full trisomy. Decisions:

1. **Runs are named and configured** (`config/runs/<name>.R`, `scripts/lib/run.R`).
   Outputs go to `results/runs/<name>/`. Baseline (Hunter-style, no
   covariates) is regenerated through the mechanism and audited against the
   archived 2026-09-04 outputs (`scripts/audit_baseline_drift.R`).
2. **Adjusted data reaches the eQTL stage as one artifact.** Script 01
   writes log2-CPM (chr21-excluded library size) with the run's covariate
   effects removed and karyotype kept; script 03 reads it. The alternative, covariates inside every stage's own
   model, was rejected as harder to audit.
3. **The adjusted run's covariates are age, sex, BMI, sample source and 19
   CyTOF cell fractions** (20 minus the reference cluster, classical
   monocytes and M-MDSCs), used directly rather than as principal
   components so the adjustment can be attributed to named cell types
   (`S1_attribution.csv`). Age, sex and sample source are the covariate set
   the HTP group's own whole-blood analyses use (Waugh 2023; Galbraith 2023
   used surrogate variables). BMI and composition are plausibly downstream
   of trisomy, so the adjusted run estimates the deviation net of those
   pathways and is reported beside, not instead of, the baseline.
4. **Mosaic T21 subjects are excluded from the adjusted run** by a cohort
   rule; translocation and unspecified (DS_T21) subjects are kept and
   flagged. Ploidy 1.5 does not hold for mosaics, and S1 panel C shows them
   apart genome-wide (9 chr21 and 58 non-chr21 genes at 5% FDR against
   full trisomy, n = 9).
5. **One intentional change to baseline:** the within-T21 eQTL fits moved
   from log2(raw count + 1) to the log2-CPM artifact so one expression scale
   runs through the pipeline. The drift audit shows every DESeq2 table and
   lane classification identical, and one eQTL call flips: CYYR1, from
   no_cis_eqtl (q 0.142) to cis_eqtl (q 0.041), because library-size
   normalization sharpens its best variant (min p 2.4e-3 to 5.9e-4). The
   baseline headline is therefore 8 DE_high cis-eQTL genes rather than 7,
   and CYYR1 leaves the set of gene-specific deviations without a detected
   cis-eQTL. This flip is recorded here for the user's review.
6. **Script 10 (log-CPM composition adjustment) is retired** to
   `scripts/archive/`; `10_compare_runs.R` compares the two full DESeq2
   runs instead.

A drift-audit lesson: the first regenerated baseline shifted four labels
of the since-retired neighborhood check because its expressed-gene pool had
been computed over all 399 count-matrix columns instead of the 397 cohort
samples. The pool (`run$partner_pool(counts, lab_ids)`, now used only by
S1) is cohort-restricted.


### Neighborhood-check null: correlation-matched, not independent (retired 2026-09-10)

Kept for the record of the retired check. The null had to be matched, not independent.
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

### eQTL controls v2: matched positive controls, several decoy sets per gene (2026-09-30)

Both standalone controls of the gene-level test were rebuilt after a
review of what the first versions showed.

**Positive controls.** The first version took the 10 chr21 genes with the
smallest GTEx whole-blood p among expressed, non-deviating genes. Those
were the largest eQTLs on the chromosome (aFC up to 5.8), nine of them in
the distal 3 Mb, several sharing a locus (GATD3A and PWP2 on one variant;
DIP2A and S100B 800 bp apart), two inside a GRCh38 false duplication
(GATD3A; KCNE1 is a deviating gene in the same list). All 10 hit the
permutation floor, which showed that the test finds very strong eQTLs and
nothing about its power for effects of the size the deviating genes carry
(GTEx aFC 0.3 to 1.4). Script 02 now picks one control per tested
deviating gene (`match_positive_controls`): among GTEx eGenes (q < 0.05,
`positive_egene_qval`) that are expressed, non-repeat and non-deviating,
the nearest in standardised (log2 |aFC|, log10 baseMean), greedily by
decreasing target |aFC|, with TSS at least 1 Mb from every other control
(`positive_min_separation`) and 100 kb from any deviating gene
(`positive_dev_separation`; a 1 Mb exclusion cut the pool from 103 to 40
and left poorer matches), never a gene in
`data/grch38_false_duplication_genes_chr21.csv`. Adjusted run: 7 of 10
matched controls detected, the same rate as the deviating genes (7 of 10).
The three missed (MRPL39, NRIP1, AF165147.1) have GTEx |aFC| 0.25 to 0.50
and, for AF165147.1, baseMean 38: the test loses power below about 0.3
log2 aFC and at low expression. AP000692.1 (CBR3's match) is detected
on a single retained variant with the opposite sign to GTEx; a weak
lncRNA eGene, kept as the match but not to be leaned on.

**Negative controls.** The first version paired each deviating gene with
one decoy set; 10 tests reused 5 sets, two landed on MAF 0.02 variants
(useless intervals), and 4 of 10 had permutation p < 0.08, which 10
correlated tests cannot interpret. Script 03 now tests each deviating
gene against `n_decoy_sets` (5) sets from all tested genes' cis variants
(deviating genes and positive controls as donors) at least 5 Mb away,
closest variant count first, common variants only (`decoy_min_maf`, GTEx
MAF 0.05), and writes one row per (gene, set) with `decoy_rank`; scripts
11 and 14 draw the rank-1 set. Adjusted run: 0 of 50 detected at FDR
0.05; 6 of 50 (12%) at nominal permutation p < 0.05 against 5% expected
(binomial p = 0.04). Three of the six are ABCC13 against decoy sets at
36, 45.9 and 46.1 Mb, so ABCC13 expression carries some chr21-wide
genotype structure; the model has no ancestry covariates (only chr21 is
genotyped, so genotype PCs would have to come from the WGS producers).
The nominal rate is a mild excess, not evidence that the FDR calls are
inflated.

**What the negatives calibrate.** A decoy best variant reaches |z| about
2.3 (median) and up to 3.4 by selection alone. Effect sizes at the best
variant are therefore not comparable across sets without the permutation
q; report the q for detection and effect sizes at the GTEx lead variant
(script 14's sensitivity table) where T21 has no winner's curse.

### Power by spike-in and the base rate on Expected-dosage genes (2026-09-30)

Two more controls added on the same branch, after the matched positive
controls (above) were tightened to strong eGenes (`positive_egene_qval`
1e-4, `positive_min_variants` 10 at the pval cut): the first matching let
in AP000692.1, a lncRNA eGene at GTEx q 0.03 with one retained variant,
"detected" in T21 with the opposite sign. Adjusted run after tightening:
8 of 10 matched positives detected (MRPL39 and NRIP1 missed; GTEx aFC
0.29 and 0.25).

**Spike-in power (script 15).** Real matched genes cannot share a gene's
own variant set, LD, allele frequencies or noise, so the power question
for a gene that was "tested, not detected" is answered by spiking its
GTEx aFC into its own permuted expression at its GTEx lead variant, on
its own T21 genotypes, and scoring the best-variant search against its
own permutation null. Adjusted run, 200 simulations per point: power at
the GTEx aFC is 1.00 for BACE2 and OLIG2 and 0.81 for ABCC13 (lead
variant MAF 0.05, so ABCC13 needs aFC 1.4 for 80%); 0.80 to 1.00 for the
seven detected genes. The aFC for 80% power is 0.27 to 0.45 for every
gene but ABCC13. So BACE2 and OLIG2 are not missed for want of power:
their GTEx eQTLs, if acting in T21 at GTEx strength, would have been
found. The spike uses the observed expression permuted, which keeps any
real eQTL variance in the noise, so the power is if anything understated.

**Base rate (script 16).** The identical test on every Expected-dosage
gene with GTEx cis variants at the cut (99 of 129 testable after the
expression join; 500 permutations): 56 of 99 detected (57%), against 7 of
10 deviating genes (70%), Fisher p = 0.51; 33 of 59 (56%) among
Expected-dosage genes in the deviating genes' baseMean range. Having a
detectable cis-eQTL is the norm for an expressed chr21 gene in this
cohort and does not distinguish the deviating genes. This is the control
the manuscript's "7 of 10" needs beside it.

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
  `chr21_lane_flow.csv` keeps them.
- **Two T21 expression-only subjects**: 304 T21 in the expression cohort,
  302 in the genotype cohort. The 2 missing-genotype subjects are
  expression-only and are silently dropped from the within-T21 regression
  step. Phrase paper text accordingly ("302 of 304").
- **lfcSE column**: chr21 combined output (`deseq2_chr21_genes_both_analyses.csv`)
  does not carry `lfcSE`; the all-genes ploidy-normalized output does.
  The archived forest-plot script joins lfcSE in from the all-genes table.
