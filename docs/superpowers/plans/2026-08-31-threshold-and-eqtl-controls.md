# Hunter Compliance, Deviation Threshold, and eQTL Controls Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bring the DE pipeline into compliance with Hunter et al.'s five prerequisites for ploidy-corrected DE, replace the arbitrary cohort-SD deviation threshold with a statistically determined chr21-internal outlier test, give the eQTL "explained" call real positive and negative controls, and add a Table 1 for the RNA-seq + WGS analysis cohort.

**Architecture:** All analysis is integrated into the **existing 00-07 pipeline scripts** — no new numbered stages. Shared statistical machinery lives in three pure-function modules under `scripts/lib/` so it can be unit-tested with testthat while the pipeline scripts stay thin orchestrators. Every script with replaced code gets a `CHANGELOG` block appended at its end. Permutation-heavy work runs locally at reduced scale by default and offloads to Alpine for the full run.

**Tech Stack:** R >= 4.2, data.table, tidyverse, DESeq2, testthat, ggplot2/patchwork. **No new package dependencies** — `mixtools` and `qvalue` are NOT installed and must not be used; BH-FDR via base `p.adjust` suffices. Alpine (CU Boulder CURC) SLURM for the full permutation run.

**Spec:** `docs/METHODS_SPEC_threshold_and_eqtl_controls.md`

**Execution order matters.** Task 1 defines the analysis cohort that Task 2 runs DE on, and Task 2 changes every ploidy-corrected number in the pipeline. Do not tune thresholds or interpret controls before both have landed.

## Global Constraints

- **No emojis** anywhere — code, comments, commit messages, output.
- **No nested loops** over data. Permutation code must use matrix algebra.
- R only for analysis code. Never commit to `main`; branch, then PR.
- **Integrate into existing pipeline scripts.** No new numbered stages. New code goes in `scripts/00`-`scripts/07` or `scripts/lib/`.
- **Every script with replaced code gets a CHANGELOG block appended at its very end**, in this form:

```r
# =============================================================================
# CHANGELOG
# =============================================================================
# YYYY-MM-DD  <what was replaced> -> <what replaced it>.
#             Reason: <why the old code was wrong>.
#             Spec: docs/METHODS_SPEC_threshold_and_eqtl_controls.md
```

- **Documentation scope:** update `CLAUDE.md` **only** where it describes the pipeline — script catalog, methodology, constants, inputs/outputs. **Do NOT update headline result numbers or figure references in `CLAUDE.md` or `README.md`.** Append decision records to `docs/REPO_STATE.md` (gitignored working notes).
- Constants that must not change value: `ALPHA = 0.01`, `ALPHA_REPRO = 0.05`, `LOW_EXPR_QUANT = 0.20`, `GTEX_PVAL_KEEP = 1e-4`, `RESTRICT_TO_PROTEIN_CODING = TRUE`.
- New constants: `OUTLIER_FDR = 0.10`, `N_PERM = 1000` (full) / `100` (local), `FDR_GENE = 0.05`. Optional: `LOW_EXPR_ABS = 30` if the low-expression label is switched from the q20 quantile to Hunter's literal absolute cutoff.
- **Never convert a label into a filter.** `low_expr` and `high_repeat` mark genes; they must not remove rows from any output table.
- **Values measured BEFORE Task 2** — expect them to move once the size-factor bug is fixed, and re-measure rather than assuming: chr21 corrected median `-0.0268`, MAD `0.2364` (all 160); after filters n = 119, median `-0.0225`, MAD `0.2257`; FDR<0.10 gives 5 genes at k `2.91`. The analysis cohort (**302 T21 + 95 Control = 397**) is a cohort fact and will not move.
- Every stage writes `*_session_info.txt` and **verifies its output files exist after writing** — a failed graphics device is otherwise silent.
- `device = "pdf"` in ggsave, never `cairo_pdf` (no X11/cairo here).

---

## Hunter et al. compliance audit

The five prerequisites from Hunter et al. (2023), audited against the current pipeline:

| # | Hunter step | Status | Where |
|---|---|---|---|
| 1 | Minimum read coverage filter (they used 30) | **IMPLEMENTED as a label** | `LOW_EXPR_QUANT = 0.20` gives a q20 baseMean cut of 25.10, flagging 32 of 160 chr21 protein-coding genes into the `Low_expression` lane. Hunter's absolute equivalent (baseMean < 30) flags 34; the 32 are a strict subset, so the rules differ by 2 genes. Hunter's *count* form (total counts >= 30) is already satisfied by every gene - the minimum is 437, median 195,547. Task 2 optionally tightens the label to the literal 30; it does NOT filter. |
| 2 | Mask repeat regions / remove multi-mapping reads at counting | **IMPLEMENTED as a label, at gene rather than read level** | `scripts/archive/process_blacklist.R` overlaps chr21 genes against ENCODE `hg38-blacklist.v2.bed.gz` regions -> `data/processed/blacklisted_genes.csv` (10 genes), combined with 9 `KNOWN_REPEAT_GENES` into the `high_repeat` flag and the `High_repeats` lane. Read-level masking is not available to us because counts arrive precomputed from Synapse. State that as a limitation; no code change. |
| 3 | Remove chr21 genes for size factor calculation | **RAW arm yes, NORM arm NO** | `01:175` does it correctly for the raw arm. `01:210` sets `normalizationFactors(dds_norm) <- norm_matrix` for the ploidy arm, which **replaces size factors entirely** — so there are none to exclude chr21 from, and no depth normalization at all. Task 2 fixes it. |
| 4 | For noisy samples, increase depth or replication | **N/A, but related** | n = 399 is large. The relevant residue is DESeq2 nulling p-values for count-outlier genes (MX1). Task 2 adds the Cook's sensitivity arm. |
| 5 | Adjust the null / normalize counts by ploidy | **YES** | The 1.5/1.0 ploidy matrix is built at `01:134-144` and applied at `01:210`. Correct in intent; Task 2 fixes how it is applied. |

**Design principle: label, do not filter.** Hunter et al. *remove* genes failing steps 1 and 2. This pipeline deliberately *labels* them instead - `low_expr` and `high_repeat` flags feeding the `Low_expression` and `High_repeats` lanes - so that genes below threshold stay in the table and remain available for later eQTL-support analysis. Every task in this plan preserves that: nothing is dropped from `chr21_lane_assignments.csv`. Where a filter is needed (estimating the robust null in Tasks 3-4), it restricts the *estimation set*, not the *reported set*, and the excluded genes carry `q_outlier = NA` rather than disappearing.

**Where we deliberately go beyond Hunter et al.:** the MAD/FDR outlier machinery (Tasks 3-4) and the eQTL controls (Tasks 6-9) have no counterpart in the paper. Hunter et al. worked with a family of four, where small `n` produced large `lfcSE` and only large effects survived `padj < 0.01` — an effect-size filter for free. At n = 399 that protection is gone, so an explicit one is required. This is an addition to their method, not a departure from it.

---

## File Structure

**Created:**

| Path | Responsibility |
|---|---|
| `scripts/lib/chr21_threshold.R` | Robust null, robust z, BH-FDR outlier test, k sensitivity |
| `scripts/lib/eqtl_fit.R` | Vectorized per-variant regression, permutation null, gene-level p, expression lookup |
| `scripts/lib/table1.R` | Continuous/categorical summary rows, two-group tests |
| `tests/testthat.R` + `tests/testthat/test-{chr21-threshold,eqtl-fit,table1}.R` | Unit tests |
| `scripts/alpine/ralpine` | Version-controlled Alpine access boundary (re-homed from fm-pdo-evaluator) |
| `scripts/alpine/permutation_controls.sbatch` | Full-scale permutation run |

**Modified (each gains a CHANGELOG block):** `scripts/01`-`04`, `scripts/07`, `scripts/download_vcf_array.sh`, `CLAUDE.md` (pipeline description only), `docs/REPO_STATE.md` (decision log).

**Deleted:** `scripts/05_chr21_distribution_panel.R` (superseded duplicate of `06_`).

---

## Task 1: Define the analysis cohort and produce Table 1

**Files:**
- Create: `scripts/lib/cohort.R`, `scripts/lib/table1.R`
- Create: `tests/testthat.R`, `tests/testthat/test-cohort.R`, `tests/testthat/test-table1.R`
- Modify: `scripts/00_preprocess_data.R` (add the WGS flag, write the cohort roster, emit Table 1; CHANGELOG at end)

**Interfaces:**
- Consumes: `data/chr21_ds_PASS.csv` and `data/chr21_ctrl_PASS.csv` (**headers only** — never read the 6 GB body here), `data/P4C_metadata_021921_Costello.txt`, optionally `data/P4C_Comorbidity_020921.tsv`.
- Produces:
  - `wgs_subjects(vcf_paths)` -> character vector of subject IDs carrying genotypes
  - `subject_id_from_labid(labid)` -> character vector
  - `analysis_cohort(meta)` -> the metadata subset defining the cohort
  - `summarize_continuous(x, digits)` -> `"median [q25-q75]"`; `summarize_categorical(x, level)` -> `"n (pct%)"`; `compare_groups(x, group)` -> `list(p, test)`
  - `data/processed/sample_metadata.csv` gains `subject_id` and `has_wgs` columns
  - `data/processed/analysis_cohort.csv` — the roster scripts 01-04 must subset to
  - `results/tables/table1_analysis_cohort.csv` / `.md`

**Cohort definition — read before writing code.**

Genotypes are used **only** for the within-T21 dosage regressions. Controls
never enter those, so requiring WGS of a control would discard 89 of 95 (94%)
for no analytic gain. The cohort is therefore asymmetric:

| arm | requirement | n |
|---|---|---|
| T21 (DS) | RNA-seq **and** WGS | **302** (of 304 with RNA-seq) |
| Control (D21) | RNA-seq only | **95** (all of them) |
| total | | **397** |

The 2 T21 subjects without genotypes are dropped from the DE analysis so that
the DE cohort and the eQTL cohort are the same people — the pipeline previously
ran DE on 304 and eQTL on 302, which made "302 of 304" phrasing necessary. Note
the DS PASS file carries 419 sample columns and the control file 14 (433 unique
subjects), so most genotyped subjects have no RNA-seq and are irrelevant here.

`subject_id` is derived from `LabID` by stripping the trailing visit suffix
(`sub("[A-Z][0-9]*$", "", LabID)`), matching script 02 line ~186. **Joining on
`RecordID` yields zero matches** — verified. Metadata fields are exactly
`RecordID, Sex, Karyotype, Event_name, LabID, Age_at_visit, BMI, Sample_source`;
there is **no** race or ethnicity field, so note its absence rather than
omitting it silently.

- [ ] **Step 1: Create the testthat runner**

Create `tests/testthat.R`:

```r
library(testthat)
test_dir("tests/testthat")
```

- [ ] **Step 2: Write the failing cohort tests**

Create `tests/testthat/test-cohort.R`:

```r
source(file.path("..", "..", "scripts", "lib", "cohort.R"))

test_that("subject_id_from_labid strips the trailing visit suffix", {
  expect_equal(subject_id_from_labid(c("HTP0001A", "HTP0002B2", "HTP0003")),
               c("HTP0001", "HTP0002", "HTP0003"))
})

test_that("wgs_subjects reads headers only and drops the fixed VCF columns", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,HTP0001A,HTP0002B",
               "chr21,100,.,A,G,.,PASS,.,GT,0/1,1/1"), f)
  expect_setequal(wgs_subjects(f), c("HTP0001", "HTP0002"))
})

test_that("wgs_subjects unions across files and de-duplicates", {
  f1 <- tempfile(fileext = ".csv"); f2 <- tempfile(fileext = ".csv")
  writeLines("CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,HTP0001A", f1)
  writeLines("CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,HTP0001B,HTP0009A", f2)
  expect_setequal(wgs_subjects(c(f1, f2)), c("HTP0001", "HTP0009"))
})

test_that("analysis_cohort keeps genotyped T21 and ALL controls", {
  meta <- data.table::data.table(
    LabID     = c("A1A", "A2A", "A3A", "C1A", "C2A"),
    Karyotype = c("T21", "T21", "T21", "Control", "Control"),
    has_wgs   = c(TRUE, TRUE, FALSE, TRUE, FALSE))
  got <- analysis_cohort(meta)
  expect_equal(sum(got$Karyotype == "T21"), 2L)      # A3A dropped, no WGS
  expect_equal(sum(got$Karyotype == "Control"), 2L)  # both kept regardless
})

test_that("analysis_cohort refuses metadata missing has_wgs", {
  meta <- data.table::data.table(LabID = "A1A", Karyotype = "T21")
  expect_error(analysis_cohort(meta), "has_wgs")
})
```

- [ ] **Step 3: Write the failing Table 1 tests**

Create `tests/testthat/test-table1.R`:

```r
source(file.path("..", "..", "scripts", "lib", "table1.R"))

test_that("summarize_continuous formats median and IQR", {
  expect_equal(summarize_continuous(c(1,2,3,4,5,6,7,8,9,10), digits = 1), "5.5 [3.2-7.8]")
})

test_that("summarize_continuous excludes missing values", {
  expect_match(summarize_continuous(c(1, 2, 3, NA, 5), digits = 1), "^3\\.0 ")
})

test_that("summarize_continuous handles an all-missing vector", {
  expect_equal(summarize_continuous(c(NA_real_, NA_real_)), "-")
})

test_that("summarize_categorical gives count and percent of non-missing", {
  x <- c("Male", "Female", "Female", "Female", NA)
  expect_equal(summarize_categorical(x, "Female"), "3 (75.0%)")
  expect_equal(summarize_categorical(x, "Male"), "1 (25.0%)")
})

test_that("summarize_categorical handles an absent level", {
  expect_equal(summarize_categorical(c("A", "A"), "B"), "0 (0.0%)")
})

test_that("compare_groups uses Wilcoxon for continuous and Fisher for categorical", {
  set.seed(20)
  cont <- compare_groups(c(rnorm(30), rnorm(30, 3)), rep(c("a", "b"), each = 30))
  expect_equal(cont$test, "Wilcoxon rank-sum")
  expect_lt(cont$p, 0.001)
  cat2 <- compare_groups(c(rep("M", 25), rep("F", 5), rep("M", 5), rep("F", 25)),
                         rep(c("a", "b"), each = 30))
  expect_equal(cat2$test, "Fisher exact")
  expect_lt(cat2$p, 0.001)
})

test_that("compare_groups refuses anything but two groups", {
  expect_error(compare_groups(1:10, rep("a", 10)), "exactly two groups")
})
```

- [ ] **Step 4: Run both test files to verify they fail**

Run: `Rscript tests/testthat.R`
Expected: FAIL — `could not find function "subject_id_from_labid"` / `"summarize_continuous"`.

- [ ] **Step 5: Write the cohort library**

Create `scripts/lib/cohort.R`:

```r
# cohort.R
#
# Analysis-cohort definition.
#
# The cohort is deliberately ASYMMETRIC. Genotypes are used only for the
# within-T21 dosage regressions, which controls never enter, so requiring WGS of
# a control would discard 89 of 95 (94%) for no analytic gain:
#
#   T21     - needs RNA-seq AND WGS  -> 302 of 304
#   Control - needs RNA-seq only     -> 95 of 95
#
# Dropping the 2 ungenotyped T21 makes the DE cohort and the eQTL cohort the
# same people. Previously DE ran on 304 and eQTL on 302.

VCF_FIXED_COLS <- c("CHROM", "POS", "ID", "REF", "ALT",
                    "QUAL", "FILTER", "INFO", "FORMAT")

#' Strip the trailing visit suffix from a LabID to get a subject ID.
#' Matches the rule in scripts/02_filter_genotypes.R (~line 186). Joining
#' genotypes on RecordID instead gives zero matches.
subject_id_from_labid <- function(labid) sub("[A-Z][0-9]*$", "", labid)

#' Subject IDs carrying chr21 genotypes, from VCF-style CSV HEADERS ONLY.
#'
#' Reads with nrows = 0 so this never touches the 6 GB body - it must stay cheap
#' enough to run in script 00, before any genotype processing.
wgs_subjects <- function(vcf_paths) {
  cols <- unlist(lapply(vcf_paths, function(f) {
    if (!file.exists(f)) {
      warning("genotype file not found, skipping: ", f)
      return(character(0))
    }
    names(data.table::fread(f, nrows = 0))
  }), use.names = FALSE)
  unique(subject_id_from_labid(setdiff(cols, VCF_FIXED_COLS)))
}

#' Subset metadata to the analysis cohort: genotyped T21 plus ALL controls.
analysis_cohort <- function(meta) {
  if (!"has_wgs" %in% names(meta)) {
    stop("metadata needs a has_wgs column; call wgs_subjects() first")
  }
  meta[(Karyotype == "T21" & has_wgs) | Karyotype != "T21"]
}
```

- [ ] **Step 6: Write the Table 1 library**

Create `scripts/lib/table1.R`:

```r
# table1.R
#
# Summary helpers for the analysis-cohort characteristics table.
#
# Continuous variables are median [IQR] with a rank-based test, not mean (SD)
# with a t-test: Age_at_visit and BMI are right-skewed. Categorical comparisons
# use Fisher exact rather than chi-square so the helpers stay valid on small
# strata.

summarize_continuous <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("-")
  q <- quantile(x, c(0.25, 0.5, 0.75), names = FALSE, type = 7)
  sprintf("%.*f [%.*f-%.*f]", digits, q[2], digits, q[1], digits, q[3])
}

summarize_categorical <- function(x, level) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("-")
  n <- sum(as.character(x) == as.character(level))
  sprintf("%d (%.1f%%)", n, 100 * n / length(x))
}

compare_groups <- function(x, group) {
  g <- factor(group)
  if (nlevels(g) != 2) stop("compare_groups needs exactly two groups; got ", nlevels(g))
  keep <- !is.na(x)
  x <- x[keep]; g <- droplevels(g[keep])
  if (nlevels(g) != 2) return(list(p = NA_real_, test = "not comparable"))
  if (is.numeric(x)) {
    list(p = suppressWarnings(wilcox.test(x ~ g)$p.value), test = "Wilcoxon rank-sum")
  } else {
    tab <- table(as.character(x), g)
    list(p = tryCatch(fisher.test(tab, simulate.p.value = nrow(tab) > 2)$p.value,
                      error = function(e) NA_real_),
         test = "Fisher exact")
  }
}
```

- [ ] **Step 7: Run the tests to verify they pass**

Run: `Rscript tests/testthat.R`
Expected: PASS.

- [ ] **Step 8: Wire the cohort definition into script 00**

In `scripts/00_preprocess_data.R`, after the metadata is matched to the count
matrix and before `sample_metadata.csv` is written:

```r
source("scripts/lib/cohort.R")

# WGS roster from the PASS file headers only - cheap, and needed here because
# script 01 must run DE on the analysis cohort, before script 02 processes any
# genotypes.
wgs <- wgs_subjects(c("data/chr21_ds_PASS.csv", "data/chr21_ctrl_PASS.csv"))
cat(sprintf("  Subjects with chr21 WGS (from headers): %d
", length(wgs)))

metadata[, subject_id := subject_id_from_labid(LabID)]
metadata[, has_wgs    := subject_id %in% wgs]
print(metadata[, .(rnaseq = .N, with_wgs = sum(has_wgs)), by = Karyotype])

cohort <- analysis_cohort(metadata)
cat(sprintf("  Analysis cohort: T21 %d, Control %d, total %d
",
            sum(cohort$Karyotype == "T21"), sum(cohort$Karyotype == "Control"),
            nrow(cohort)))
stopifnot(sum(cohort$Karyotype == "T21") == 302,
          sum(cohort$Karyotype == "Control") == 95)
fwrite(cohort, "data/processed/analysis_cohort.csv")
```

Keep writing the **full** `sample_metadata.csv` (all 399 rows, now carrying
`subject_id` and `has_wgs`) — downstream code that wants the cohort reads
`analysis_cohort.csv`. Do not silently shrink the metadata file.

- [ ] **Step 9: Add Table 1 to script 00**

After the block above:

```r
source("scripts/lib/table1.R")

t21a <- cohort[Karyotype == "T21"]; ctla <- cohort[Karyotype == "Control"]

row_cont <- function(label, var, digits = 1) {
  cmp <- compare_groups(cohort[[var]], cohort$Karyotype)
  data.table(characteristic = label, level = "median [IQR]",
             t21 = summarize_continuous(t21a[[var]], digits),
             control = summarize_continuous(ctla[[var]], digits),
             overall = summarize_continuous(cohort[[var]], digits),
             p_value = cmp$p, test = cmp$test,
             n_missing = sum(is.na(cohort[[var]])))
}
row_cat <- function(label, var) {
  cmp  <- compare_groups(cohort[[var]], cohort$Karyotype)
  levs <- sort(unique(as.character(cohort[[var]][!is.na(cohort[[var]])])))
  rbindlist(lapply(seq_along(levs), function(i)
    data.table(characteristic = if (i == 1) label else "", level = levs[i],
               t21 = summarize_categorical(t21a[[var]], levs[i]),
               control = summarize_categorical(ctla[[var]], levs[i]),
               overall = summarize_categorical(cohort[[var]], levs[i]),
               p_value = if (i == 1) cmp$p else NA_real_,
               test = if (i == 1) cmp$test else "",
               n_missing = if (i == 1) sum(is.na(cohort[[var]])) else NA_integer_)))
}

tbl1 <- rbindlist(list(
  data.table(characteristic = "Subjects in analysis cohort", level = "n",
             t21 = as.character(nrow(t21a)), control = as.character(nrow(ctla)),
             overall = as.character(nrow(cohort)),
             p_value = NA_real_, test = "", n_missing = NA_integer_),
  data.table(characteristic = "Excluded: RNA-seq without WGS", level = "n",
             t21 = as.character(sum(metadata$Karyotype == "T21" & !metadata$has_wgs)),
             control = "0 (not required)",
             overall = as.character(sum(metadata$Karyotype == "T21" & !metadata$has_wgs)),
             p_value = NA_real_, test = "", n_missing = NA_integer_),
  row_cont("Age at visit (years)", "Age_at_visit"),
  row_cont("BMI", "BMI"),
  row_cat("Sex", "Sex"),
  row_cat("Sample source", "Sample_source"),
  row_cat("Study visit", "Event_name")))

# Optional comorbidity block. Long format with ONE comment line before the
# header: RecordID, Condition, HasCondition, Age.group, min_Age, max_Age.
como_path <- "data/P4C_Comorbidity_020921.tsv"
if (file.exists(como_path)) {
  como <- fread(como_path, skip = 1)
  if (all(c("RecordID", "Condition", "HasCondition") %in% names(como))) {
    como <- como[RecordID %in% cohort$RecordID]
    top  <- como[, .(n_with = sum(HasCondition == 1, na.rm = TRUE)), by = Condition][
      order(-n_with)][seq_len(min(5, .N))]
    wide <- dcast(como, RecordID ~ Condition, value.var = "HasCondition")
    ann  <- merge(cohort[, .(RecordID, Karyotype)], wide, by = "RecordID")
    tbl1 <- rbindlist(list(tbl1, rbindlist(lapply(top$Condition, function(cn) {
      cmp <- compare_groups(as.character(ann[[cn]]), ann$Karyotype)
      data.table(characteristic = cn, level = "n (%) with condition",
                 t21 = summarize_categorical(ann[Karyotype == "T21"][[cn]], 1),
                 control = summarize_categorical(ann[Karyotype == "Control"][[cn]], 1),
                 overall = summarize_categorical(ann[[cn]], 1),
                 p_value = cmp$p, test = cmp$test, n_missing = sum(is.na(ann[[cn]])))
    }))))
  } else {
    cat("  comorbidity file present but unexpected columns; skipping that block\n")
  }
}

tbl1[, p_value := ifelse(is.na(p_value), "", format.pval(p_value, digits = 2, eps = 1e-4))]
fwrite(tbl1, "results/tables/table1_analysis_cohort.csv")
writeLines(c(
  "# Table 1. Characteristics of the analysis cohort", "",
  sprintf("T21 subjects require both whole-blood RNA-seq and chr21 WGS (n = %d of %d).",
          nrow(t21a), sum(metadata$Karyotype == "T21")),
  sprintf("Controls require RNA-seq only (n = %d), because genotypes are used solely",
          nrow(ctla)),
  "for the within-T21 dosage regressions, which controls do not enter.",
  "Race and ethnicity are not recorded in the available metadata.", "",
  paste("|", paste(names(tbl1), collapse = " | "), "|"),
  paste("|", paste(rep("---", ncol(tbl1)), collapse = " | "), "|"),
  apply(tbl1, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))),
  "results/tables/table1_analysis_cohort.md")
for (f in c("results/tables/table1_analysis_cohort.csv",
            "results/tables/table1_analysis_cohort.md")) {
  if (!file.exists(f)) stop("failed to write ", f)
  cat("  Wrote ", f, "\n", sep = "")
}
print(tbl1)
```

- [ ] **Step 10: Run script 00 and assert the cohort**

Run: `Rscript scripts/00_preprocess_data.R 2>&1 | tail -40`
Expected: `Analysis cohort: T21 302, Control 95, total 397`. The `stopifnot` in Step 8 is the guard; if it fires, the `subject_id` derivation or the header read is wrong — fix that before trusting any row of Table 1.

- [ ] **Step 11: Sanity-check one row by hand**

```bash
Rscript -e '
suppressPackageStartupMessages(library(data.table))
a <- fread("data/processed/analysis_cohort.csv")[Karyotype == "T21"]
cat("T21 n:", nrow(a), "\n")
cat("age median [IQR]:", sprintf("%.1f [%.1f-%.1f]", median(a$Age_at_visit),
    quantile(a$Age_at_visit,.25), quantile(a$Age_at_visit,.75)), "\n")
cat("female:", sum(a$Sex == "Female"), sprintf("(%.1f%%)", 100*mean(a$Sex=="Female")), "\n")
'
```

Expected: matches the T21 column of Table 1 exactly.

- [ ] **Step 12: Append the CHANGELOG block to script 00 and commit**

```r
# =============================================================================
# CHANGELOG
# =============================================================================
# 2026-08-31  ADDED the analysis-cohort definition and Table 1.
#             The cohort is asymmetric by design: T21 require RNA-seq AND chr21
#             WGS (302 of 304), Controls require RNA-seq only (95 of 95).
#             Reason: genotypes are used only for the within-T21 dosage
#             regressions, which controls never enter, so requiring WGS of a
#             control would discard 89 of 95 for no analytic gain. Dropping the
#             2 ungenotyped T21 makes the DE cohort and the eQTL cohort the same
#             people; previously DE ran on 304 and eQTL on 302.
#             The WGS roster is read from the PASS file HEADERS only, so this
#             stays cheap enough to run before any genotype processing.
#             sample_metadata.csv keeps all 399 rows and gains subject_id and
#             has_wgs; data/processed/analysis_cohort.csv carries the 397-row
#             roster that scripts 01-04 subset to.
#             Spec: docs/METHODS_SPEC_threshold_and_eqtl_controls.md
```

```bash
git checkout -b hunter-threshold-eqtl-controls
git add scripts/lib/cohort.R scripts/lib/table1.R tests/testthat.R \
        tests/testthat/test-cohort.R tests/testthat/test-table1.R \
        scripts/00_preprocess_data.R
git commit -m "Define the analysis cohort (302 T21 with WGS, 95 controls) and add Table 1"
```

---

## Task 2: Hunter compliance — size-factor fix and Cook's arm in script 01

**This task changes every ploidy-corrected number downstream. It must land before Tasks 3-9, and after Task 1, whose cohort roster it consumes.**

**Files:**
- Modify: `scripts/01_deseq2_analysis.R` (size-factor call ~210; results extraction ~217; CHANGELOG at end)
- Optionally modify: `scripts/04_chr21_lane_assignment.R` (low-expression label threshold only)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `results/tables/deseq2_cooks_diagnostics.csv` (`EnsemblID, padj_default, padj_nocooks, norm_log2FC, baseMean, max_cooks, n_samples_over_cutoff`); `results/tables/deseq2_chr21_genes_both_analyses_nocooks.csv` matching the schema of the default chr21 table. Existing outputs keep their schemas; their **values change**.

- [ ] **Step 1: Align the low-expression LABEL with Hunter's threshold (optional, no filtering)**

Hunter et al. step 1 is already satisfied. Their count form (total counts >= 30)
excludes nothing here - the minimum on chr21 protein-coding is 437 total counts,
median 195,547. Their baseMean form flags 34 of 160 genes, against 32 under the
current q20 quantile rule, and the 32 are a strict subset. So this step is a
2-gene refinement for literal compliance, **not** a filter.

**Do not add a pre-DESeq2 count filter.** Removing genes would empty the
`Low_expression` and `High_repeats` lanes that the Sankey depends on and would
discard genes needed for later eQTL-support analysis.

If literal compliance is wanted, in `scripts/04_chr21_lane_assignment.R` replace
the quantile threshold with an absolute one, keeping it a label:

```r
# Hunter et al. step 1, absolute form. Kept as a LABEL, not a filter: genes
# below it stay in the table and flow to the Low_expression lane, so they remain
# available for eQTL-support analysis. The previous rule was the q20 baseMean
# quantile (25.10 here), which flags 32 genes; this flags 34, a strict superset.
LOW_EXPR_ABS <- 30
basemean_threshold <- LOW_EXPR_ABS
```

Otherwise leave the quantile rule as-is and record in the CHANGELOG that step 1
was audited and found already satisfied. Either way, make no change to script 01.

- [ ] **Step 1b: Subset the DE analysis to the analysis cohort**

Script 01 currently runs DE on all 399 samples. It must run on the 397-row
cohort from Task 1 so the DE and eQTL analyses cover the same people. Insert
immediately after `count_matrix` and `metadata` are loaded, before the ploidy
matrix is built:

```r
# Analysis cohort from Task 1: T21 need RNA-seq AND WGS (302 of 304); Controls
# need RNA-seq only (95 of 95), because genotypes are used solely for the
# within-T21 dosage regressions that controls never enter. Running DE on the
# same people as the eQTL step removes the old "302 of 304" mismatch.
cohort <- read_csv("data/processed/analysis_cohort.csv", show_col_types = FALSE)
keep_samples <- colnames(count_matrix) %in% cohort$LabID
cat(sprintf("  Analysis cohort: %d of %d samples (T21 %d, Control %d)\n",
            sum(keep_samples), ncol(count_matrix),
            sum(cohort$Karyotype == "T21"), sum(cohort$Karyotype == "Control")))
stopifnot(sum(keep_samples) == 397)
count_matrix <- count_matrix[, keep_samples, drop = FALSE]
metadata     <- metadata[metadata$LabID %in% cohort$LabID, ]
```

The `norm_matrix` construction and its `stopifnot` dimension checks come after
this, so they will catch any mis-ordering.

- [ ] **Step 2: Fix the size-factor calculation in the ploidy arm**

Replace `scripts/01_deseq2_analysis.R:210`:

```r
normalizationFactors(dds_norm) <- norm_matrix
```

with:

```r
# Hunter et al. steps 3 and 5 together. estimateSizeFactors(normMatrix = ...)
# computes library-size factors GIVEN the ploidy matrix and folds both into the
# normalization factors, and controlGenes restricts that estimation to non-chr21
# genes.
#
# The previous code assigned normalizationFactors(dds_norm) <- norm_matrix
# directly. normalizationFactors are the COMPLETE per-gene-per-sample divisors
# and replace size factors outright, and DESeq() skips estimateSizeFactors()
# when they are already set - so that arm ran with NO library-size normalization
# at all (norm_matrix holds only 1.0 and 1.5). Library sizes here span 21.1M to
# 52.3M reads, 2.47x max/min, CV 12.2%.
dds_norm <- estimateSizeFactors(dds_norm,
                                normMatrix   = norm_matrix,
                                controlGenes = non_chr21_genes_norm)
```

Define `non_chr21_genes_norm` just above it, against `dds_norm`'s own rownames:

```r
non_chr21_genes_norm <- !rownames(dds_norm) %in% chr21_genes
```

- [ ] **Step 3: Assert the fix took effect**

Immediately after the `estimateSizeFactors` call:

```r
# Guard: normalization factors must now vary across samples. If they only ever
# take the values in norm_matrix, size factors were not estimated and the arm is
# back to having no depth correction.
nf <- normalizationFactors(dds_norm)
stopifnot(!all(nf %in% c(1.0, 1.5)))
cat(sprintf("  Normalization factors: range %.3f-%.3f, %d distinct column medians\n",
            min(nf), max(nf), length(unique(round(apply(nf, 2, median), 6)))))
```

- [ ] **Step 4: Add the Cook's-distance sensitivity arm**

Immediately after `results_norm <- results(dds_norm, name = "karyotype_T21_vs_Control")`:

```r
# Hunter et al. step 4 residue. DESeq2 nulls the p-value (not just padj) for
# genes with an extreme count outlier. MX1 is affected: baseMean 9017,
# norm_log2FC 0.819, stat 5.63 - one of the largest deviations on the
# chromosome - yet pvalue NA, so the pipeline's !is.na(norm_padj) filter drops
# it. Emit both arms so the effect is visible.
results_norm_nocooks <- results(dds_norm, name = "karyotype_T21_vs_Control",
                                cooksCutoff = FALSE)

cooks_mat <- assays(dds_norm)[["cooks"]]
m_par     <- ncol(attr(dds_norm, "modelMatrix"))
cooks_cut <- qf(0.99, m_par, ncol(dds_norm) - m_par)
cooks_diag <- data.frame(
  EnsemblID             = rownames(results_norm),
  padj_default          = results_norm$padj,
  padj_nocooks          = results_norm_nocooks$padj,
  norm_log2FC           = results_norm$log2FoldChange,
  baseMean              = results_norm$baseMean,
  max_cooks             = apply(cooks_mat, 1, max, na.rm = TRUE),
  n_samples_over_cutoff = rowSums(cooks_mat > cooks_cut, na.rm = TRUE),
  stringsAsFactors      = FALSE)
cat(sprintf("  Cook's filtering nulled %d genes; cooksCutoff=FALSE recovers them\n",
            sum(is.na(cooks_diag$padj_default) & !is.na(cooks_diag$padj_nocooks))))
write_csv(cooks_diag, "results/tables/deseq2_cooks_diagnostics.csv")
cat("  Saved: results/tables/deseq2_cooks_diagnostics.csv\n")
```

And after `write_csv(chr21_combined, "results/tables/deseq2_chr21_genes_both_analyses.csv")`:

```r
# Parallel chr21 table from the cooksCutoff=FALSE arm, for the MX1 sensitivity
# analysis. Downstream scripts read the default arm.
chr21_nocooks <- chr21_combined
idx <- match(chr21_nocooks$EnsemblID, rownames(results_norm_nocooks))
chr21_nocooks$norm_pvalue <- results_norm_nocooks$pvalue[idx]
chr21_nocooks$norm_padj   <- results_norm_nocooks$padj[idx]
write_csv(chr21_nocooks, "results/tables/deseq2_chr21_genes_both_analyses_nocooks.csv")
cat("  Saved: results/tables/deseq2_chr21_genes_both_analyses_nocooks.csv\n")
```

- [ ] **Step 5: Append the CHANGELOG block**

```r
# =============================================================================
# CHANGELOG
# =============================================================================
# 2026-08-31  Hunter et al. (2023) compliance pass, plus cohort restriction.
#
#  [cohort]   RESTRICTED DE from all 399 samples to the 397-row analysis cohort
#             (302 T21 with WGS + 95 Controls, from
#             data/processed/analysis_cohort.csv). Reason: DE previously ran on
#             304 T21 while the eQTL step ran on 302, forcing "302 of 304"
#             phrasing; the two now cover the same people. Controls are NOT
#             required to have WGS - genotypes are used only for the within-T21
#             dosage regressions.
#
#  [step 3+5] REPLACED `normalizationFactors(dds_norm) <- norm_matrix` with
#             `estimateSizeFactors(dds_norm, normMatrix = norm_matrix,
#              controlGenes = non_chr21_genes_norm)`.
#             Reason: normalizationFactors are the complete per-gene-per-sample
#             divisors and replace size factors outright; DESeq() skips
#             estimateSizeFactors() when they are set. Since norm_matrix held
#             only 1.0 and 1.5, the ploidy arm ran with NO library-size
#             normalization, on samples spanning 21.1M-52.3M reads (2.47x
#             max/min, CV 12.2%). The raw arm was already correct (line ~175).
#             This changes every norm_log2FC and norm_padj in the pipeline.
#
#  [step 1]   AUDITED, no filter added. Already satisfied: the absolute count
#             form (total >= 30) excludes nothing (min 437, median 195,547), and
#             the q20 baseMean rule (25.10, flags 32 of 160) is a strict subset
#             of Hunter's baseMean < 30 (flags 34). Kept as a LABEL feeding the
#             Low_expression lane - genes below threshold stay in the table for
#             later eQTL-support analysis.
#
#  [step 4]   ADDED a cooksCutoff=FALSE sensitivity arm and Cook's diagnostics.
#             Reason: MX1 and 216 other genes have p-values nulled by outlier
#             filtering; the pipeline dropped them silently.
#
#  [step 2]   AUDITED, no change. Implemented as a LABEL at gene level:
#             scripts/archive/process_blacklist.R overlaps chr21 genes against
#             ENCODE hg38-blacklist.v2 regions (10 genes), combined with 9
#             KNOWN_REPEAT_GENES into the high_repeat flag and the High_repeats
#             lane. Read-level masking is not available because counts arrive
#             precomputed from Synapse. State as a manuscript limitation.
#
#             Spec: docs/METHODS_SPEC_threshold_and_eqtl_controls.md
```

- [ ] **Step 6: Run and record the impact**

```bash
cp results/tables/deseq2_chr21_genes_both_analyses.csv /tmp/chr21_before_hunter.csv
Rscript scripts/01_deseq2_analysis.R 2>&1 | tail -40
Rscript -e '
suppressPackageStartupMessages({library(data.table)})
a <- fread("/tmp/chr21_before_hunter.csv"); b <- fread("results/tables/deseq2_chr21_genes_both_analyses.csv")
m <- merge(a[, .(EnsemblID, old = norm_log2FC)], b[, .(EnsemblID, new = norm_log2FC)], by = "EnsemblID")
cat("chr21 genes compared:", nrow(m), "\n")
cat("median norm_log2FC:  old", round(median(m$old, na.rm=TRUE), 4),
    " new", round(median(m$new, na.rm=TRUE), 4), "\n")
cat("mean |change|:", round(mean(abs(m$new - m$old), na.rm=TRUE), 4),
    "  max:", round(max(abs(m$new - m$old), na.rm=TRUE), 4), "\n")
cat("correlation:", round(cor(m$old, m$new, use="complete.obs"), 4), "\n")
'
```

Record all of it. A large shift is expected and is the point; a *zero* shift means the fix did not take effect and Step 3's guard should have fired.

- [ ] **Step 7: Commit**

```bash
git add scripts/01_deseq2_analysis.R
git commit -m "Restrict DE to the analysis cohort; fix ploidy-arm size factors; add Cook's arm"
```

---

## Task 3: Robust chr21 threshold library

**Files:**
- Create: `scripts/lib/chr21_threshold.R`, `tests/testthat/test-chr21-threshold.R`, `tests/testthat.R`

**Interfaces:**
- Consumes: nothing (leaf library).
- Produces: `chr21_null(log2fc)` -> `list(center, scale, n)`; `robust_z(log2fc, null)` -> numeric; `outlier_fdr(z)` -> BH q-values; `effective_k(z, q, alpha)` -> numeric(1); `k_sensitivity(z, k_values)` -> `data.frame(k, n_genes, pct_genes, expected_pct_normal)`; `n_missing(x)` -> integer(1).

- [ ] **Step 1: Create the testthat runner**

Create `tests/testthat.R`:

```r
library(testthat)
test_dir("tests/testthat")
```

- [ ] **Step 2: Write the failing tests**

Create `tests/testthat/test-chr21-threshold.R`:

```r
source(file.path("..", "..", "scripts", "lib", "chr21_threshold.R"))

test_that("chr21_null returns median and scaled MAD", {
  set.seed(1)
  x <- rnorm(500, mean = 2, sd = 3)
  n <- chr21_null(x)
  expect_equal(n$center, median(x))
  expect_equal(n$scale, mad(x))
  expect_equal(n$n, 500L)
})

test_that("chr21_null drops NA and refuses tiny inputs", {
  expect_equal(chr21_null(c(1, 2, NA, 4, 5, 6, 7, 8, 9, 10, 11))$n, 10L)
  expect_error(chr21_null(c(1, 2, 3)), "at least 10")
})

test_that("chr21_null refuses a degenerate scale", {
  expect_error(chr21_null(rep(2, 50)), "scale")
})

test_that("robust_z centres and scales", {
  x <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  expect_equal(robust_z(x, list(center = 4, scale = 2, n = 10L)), (x - 4) / 2)
})

test_that("MAD scale resists outliers far better than SD", {
  set.seed(2)
  clean <- rnorm(160)
  dirty <- c(clean[1:155], 20, -20, 25, -25, 30)
  expect_gt(abs(sd(dirty) - sd(clean)) / sd(clean),
            5 * abs(mad(dirty) - mad(clean)) / mad(clean))
})

test_that("outlier_fdr finds injected outliers and spares a clean null", {
  set.seed(3)
  clean <- rnorm(200)
  expect_equal(sum(outlier_fdr(robust_z(clean, chr21_null(clean))) < 0.10), 0)
  spiked <- c(rnorm(195), 8, -8, 9, -9, 10)
  expect_equal(sum(outlier_fdr(robust_z(spiked, chr21_null(spiked))) < 0.10), 5)
})

test_that("effective_k is the smallest |z| that survives the FDR cut", {
  z <- c(0.5, 1, 2, 3, 4, 5)
  q <- c(0.9, 0.8, 0.5, 0.2, 0.02, 0.001)
  expect_equal(effective_k(z, q, 0.05), 4)
  expect_true(is.na(effective_k(z, q, 1e-6)))
})

test_that("k_sensitivity reports counts and the normal expectation", {
  set.seed(4)
  z <- rnorm(1000)
  s <- k_sensitivity(z, c(1, 2, 3))
  expect_equal(nrow(s), 3L)
  expect_equal(s$n_genes[1], sum(abs(z) >= 1))
  expect_equal(s$expected_pct_normal[2], 100 * 2 * pnorm(-2))
})

test_that("n_missing counts NA", {
  expect_equal(n_missing(c(1, NA, 3, NA)), 2L)
})
```

- [ ] **Step 3: Run to verify failure**

Run: `Rscript tests/testthat.R`
Expected: FAIL — `could not find function "chr21_null"`.

- [ ] **Step 4: Write the library**

Create `scripts/lib/chr21_threshold.R`:

```r
# chr21_threshold.R
#
# Robust outlier test for chr21 ploidy-corrected deviations.
#
# This machinery has no counterpart in Hunter et al. and is a deliberate
# addition, not a departure: their family-of-four design got an effect-size
# filter for free (small n -> large lfcSE -> only large effects cleared
# padj < 0.01). At n = 399 that protection is gone.
#
# The null is estimated from the chr21 genes themselves, not from non-chr21
# genes: ploidy normalization only acts on chr21, so a diploid gene's
# "corrected" log2FC is a different quantity and cannot be the reference.
#
# Scale is MAD rather than SD because the chr21 SD is inflated by the very
# outliers being detected (SD/MAD = 1.49; removing 5 of 160 genes moves SD 19%
# and MAD 1%, and those 5 include OLIG2, a primary candidate).
#
# Callers MUST apply the low-expression and high-repeat filters before
# estimating the null - log2FC variance scales with expression.

MIN_NULL_N <- 10

n_missing <- function(x) sum(is.na(x))

#' Estimate the robust null from chr21 corrected log2 fold changes.
#' @return list(center, scale, n); scale is mad(), already x1.4826.
chr21_null <- function(log2fc) {
  stopifnot(is.numeric(log2fc))
  x <- log2fc[!is.na(log2fc)]
  if (length(x) < MIN_NULL_N) {
    stop("need at least ", MIN_NULL_N, " non-missing values to estimate a null; got ",
         length(x))
  }
  s <- mad(x)
  if (!is.finite(s) || s <= 0) {
    stop("degenerate null scale (mad = ", s, "); the input has no spread")
  }
  list(center = median(x), scale = s, n = length(x))
}

robust_z <- function(log2fc, null) {
  stopifnot(is.list(null), all(c("center", "scale") %in% names(null)))
  if (!is.finite(null$scale) || null$scale <= 0) stop("null scale must be positive")
  (log2fc - null$center) / null$scale
}

#' Two-sided BH q-values. The chr21 bulk is normal once the 10% most extreme
#' are trimmed (Shapiro-Wilk p = 0.17 vs 6e-6 untrimmed), which licenses a
#' normal reference. Re-check that after Task 2 changes the values.
outlier_fdr <- function(z) p.adjust(2 * pnorm(-abs(z)), method = "BH")

effective_k <- function(z, q, alpha) {
  sel <- !is.na(q) & !is.na(z) & q < alpha
  if (!any(sel)) return(NA_real_)
  min(abs(z[sel]))
}

#' Gene counts across k, with the normal-null expectation. Report this instead
#' of a single count: the DE_low set is threshold-driven while DE_high is not.
k_sensitivity <- function(z, k_values = c(1, 1.5, 2, 2.5, 3, 3.5)) {
  z <- z[!is.na(z)]
  n_genes <- vapply(k_values, function(k) sum(abs(z) >= k), integer(1))
  data.frame(k = k_values, n_genes = n_genes,
             pct_genes = 100 * n_genes / length(z),
             expected_pct_normal = 100 * 2 * pnorm(-k_values))
}
```

- [ ] **Step 5: Run to verify pass**

Run: `Rscript tests/testthat.R`
Expected: PASS.

- [ ] **Step 6: Re-measure the null on the post-Task-1 data**

```bash
Rscript -e '
source("scripts/lib/chr21_threshold.R")
suppressPackageStartupMessages({library(readr); library(dplyr)})
n <- read_csv("results/tables/deseq2_all_genes_ploidy_normalized.csv", show_col_types = FALSE) %>%
  filter(Gene_type == "protein_coding", Chr == "chr21", !is.na(log2FoldChange))
nl <- chr21_null(n$log2FoldChange)
cat(sprintf("post-Task-1 chr21 null: center %.4f  MAD %.4f  n %d\n", nl$center, nl$scale, nl$n))
cat("pre-Task-1 reference was center -0.0268  MAD 0.2364  n 160\n")
b <- n$log2FoldChange
z <- abs(b - median(b)) / mad(b)
cat(sprintf("bulk normality after trimming top 10%%: Shapiro p = %.3g\n",
            shapiro.test(b[z < quantile(z, 0.90)])$p.value))
'
```

Record the new values. If the trimmed-bulk Shapiro p is now below 0.05, the normal reference in `outlier_fdr` is no longer justified — stop and reconsider before proceeding.

- [ ] **Step 7: Commit**

```bash
git add scripts/lib/chr21_threshold.R tests/testthat.R tests/testthat/test-chr21-threshold.R
git commit -m "Add robust chr21-internal threshold library with FDR outlier test"
```

---

## Task 4: Replace the magnitude filter in scripts 02 and 04

**Files:**
- Modify: `scripts/02_filter_genotypes.R` (constants; target-gene block; CHANGELOG)
- Modify: `scripts/04_chr21_lane_assignment.R` (constants; magnitude block; column list; CHANGELOG)

**Interfaces:**
- Consumes: everything from Task 3.
- Produces: `chr21_lane_assignments.csv` gains `dev_z`, `q_outlier`; `passes_magnitude_filter` becomes FDR-driven; `deviation_vs_cohort_sd` is dropped. New `results/tables/chr21_k_sensitivity.csv`.

- [ ] **Step 1: Replace the constant in script 02**

```r
# Deviation threshold. The old rule (|norm_log2FC| >= 1.0 * SD of non-chr21
# genes) used the wrong reference - ploidy normalization does not act on diploid
# genes - and 1 SD selects the top ~third of any distribution. The null is now
# chr21-internal median/MAD and the cut is an FDR on robust z.
OUTLIER_FDR <- 0.10
```

- [ ] **Step 2: Replace the target-gene selection in script 02**

Replace the `cohort_sd` computation and the `target_genes <- chr21[...]` block with:

```r
source("scripts/lib/chr21_threshold.R")

# Eligibility filters run BEFORE the null is estimated: log2FC variance scales
# with expression, so near-zero-count genes would otherwise set the scale.
basemean_threshold <- quantile(chr21$baseMean, LOW_EXPR_QUANT, na.rm = TRUE)
eligible <- chr21[baseMean >= basemean_threshold &
                    !(Gene_name %in% high_repeat_genes) &
                    !is.na(norm_log2FC)]
cat(sprintf("  Eligible for the null (expressed, non-repeat): %d of %d\n",
            nrow(eligible), nrow(chr21)))

null <- chr21_null(eligible$norm_log2FC)
cat(sprintf("  chr21 null: center %.4f  MAD %.4f  (n = %d)\n",
            null$center, null$scale, null$n))

eligible[, dev_z := robust_z(norm_log2FC, null)]
eligible[, q_outlier := outlier_fdr(dev_z)]
cat(sprintf("  Outlier test at FDR < %.2f: %d genes (effective k = %.2f)\n",
            OUTLIER_FDR, sum(eligible$q_outlier < OUTLIER_FDR, na.rm = TRUE),
            effective_k(eligible$dev_z, eligible$q_outlier, OUTLIER_FDR)))

# Composite rule: real deviation (padj) AND unusually large deviation (FDR).
target_genes <- eligible[!is.na(norm_padj) & norm_padj < ALPHA_DE &
                           !is.na(q_outlier) & q_outlier < OUTLIER_FDR]
```

- [ ] **Step 3: Append the CHANGELOG block to script 02**

```r
# =============================================================================
# CHANGELOG
# =============================================================================
# 2026-08-31  REPLACED the cohort-SD magnitude filter
#             (abs(norm_log2FC) >= MAGNITUDE_THRESHOLD * sd(non-chr21 log2FC),
#             threshold 1.0) with an FDR-controlled robust outlier test against
#             a chr21-internal median/MAD null (scripts/lib/chr21_threshold.R,
#             OUTLIER_FDR = 0.10).
#             Reason: ploidy normalization does not act on diploid genes
#             (mean |raw - norm| 0.0048 off chr21 vs 0.583 on it), so their
#             spread measured a different quantity; and a 1-SD cut selects the
#             top ~third of any distribution (18.3% of non-chr21 genes cleared
#             it themselves, vs 20.6% of chr21 - binomial p = 0.26). The null is
#             now estimated AFTER the expression and repeat filters, because
#             log2FC variance scales with counts.
#             Spec: docs/METHODS_SPEC_threshold_and_eqtl_controls.md
```

- [ ] **Step 4: Run script 02 and record the null**

Run: `Rscript scripts/02_filter_genotypes.R 2>&1 | head -30`
Expected: the null values printed. Compare against Task 3 Step 6 — they must match.

- [ ] **Step 5: Apply the same replacement in script 04**

Replace `MAGNITUDE_THRESHOLD <- 1.0` with `OUTLIER_FDR <- 0.10`, and replace the `non_chr21_sd` / `deviation_vs_cohort_sd` / `passes_magnitude_filter` block with:

```r
source("scripts/lib/chr21_threshold.R")

# Null from eligible genes only; z and q reported for all genes so the table
# stays complete. q is NA for ineligible genes, which fails the filter by
# construction - they are caught by the High_repeats / Low_expression lanes.
eligible_idx <- !m$low_expr & !m$high_repeat & !is.na(m$norm_log2FC)
null <- chr21_null(m$norm_log2FC[eligible_idx])
cat(sprintf("  chr21 null: center %.4f  MAD %.4f  (n = %d eligible genes)\n",
            null$center, null$scale, null$n))

m[, dev_z := robust_z(norm_log2FC, null)]
m[, q_outlier := NA_real_]
m[eligible_idx, q_outlier := outlier_fdr(dev_z)]
m[, passes_magnitude_filter := !is.na(q_outlier) & q_outlier < OUTLIER_FDR]

cat(sprintf("  Outlier test at FDR < %.2f: %d genes (effective k = %.2f)\n",
            OUTLIER_FDR, sum(m$passes_magnitude_filter),
            effective_k(m$dev_z, m$q_outlier, OUTLIER_FDR)))

sens <- k_sensitivity(m$dev_z[eligible_idx])
fwrite(sens, "results/tables/chr21_k_sensitivity.csv")
cat("  Wrote results/tables/chr21_k_sensitivity.csv\n")
print(sens)
```

- [ ] **Step 6: Update the written column list in script 04**

Replace `"deviation_vs_cohort_sd"` with `"dev_z", "q_outlier"`, keeping `deviation_magnitude` and `passes_magnitude_filter`.

- [ ] **Step 7: Append the CHANGELOG block to script 04**

Same as Step 3, with this first line:

```r
# 2026-08-31  REPLACED the cohort-SD magnitude filter with the chr21-internal
#             FDR outlier test; dropped column deviation_vs_cohort_sd in favour
#             of dev_z and q_outlier; added chr21_k_sensitivity.csv.
```

- [ ] **Step 8: Run 03 and 04, inspect the new gene set**

```bash
Rscript scripts/03_t21_dosage_boxplots.R && Rscript scripts/04_chr21_lane_assignment.R
Rscript -e '
suppressPackageStartupMessages(library(readr))
x <- read_csv("results/tables/chr21_lane_assignments.csv", show_col_types = FALSE)
print(table(x$sig_lane))
print(as.data.frame(subset(x, sig_lane %in% c("DE_low","DE_high"),
                           select = c(Gene_name, sig_lane, eqtl_lane, dev_z, q_outlier))))
'
```

Record the gene list. Do not compare it against pre-Task-1 expectations — those numbers came from the unnormalized arm.

- [ ] **Step 9: Commit**

```bash
git add scripts/02_filter_genotypes.R scripts/04_chr21_lane_assignment.R
git commit -m "Replace cohort-SD magnitude filter with chr21-internal FDR outlier test"
```

---

## Task 5: Vectorized eQTL regression library, and refactor script 03 onto it

**Files:**
- Create: `scripts/lib/eqtl_fit.R`, `tests/testthat/test-eqtl-fit.R`
- Modify: `scripts/03_t21_dosage_boxplots.R` (regression block ~127-160)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `fit_variants(G, e)` -> `data.frame(slope, se, t, p)`, one row per column of `n x m` matrix `G`; `perm_min_p(G, e, n_perm, seed)` -> numeric length `n_perm`; `gene_level_p(min_p_obs, min_p_perm)` -> numeric(1); `expr_of_gene(gene_name, subject_ids, counts, meta_t21)` -> numeric length `length(subject_ids)`.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-eqtl-fit.R`:

```r
source(file.path("..", "..", "scripts", "lib", "eqtl_fit.R"))

test_that("fit_variants matches lm() for a single variant", {
  set.seed(10)
  n <- 200
  g <- rbinom(n, 3, 0.3)
  e <- 0.7 * g + rnorm(n)
  got <- fit_variants(matrix(g, ncol = 1), e)
  ref <- summary(lm(e ~ g))$coefficients["g", ]
  expect_equal(got$slope[1], unname(ref["Estimate"]), tolerance = 1e-8)
  expect_equal(got$se[1],    unname(ref["Std. Error"]), tolerance = 1e-8)
  expect_equal(got$p[1],     unname(ref["Pr(>|t|)"]), tolerance = 1e-8)
})

test_that("fit_variants matches lm() for every column of a matrix", {
  set.seed(11)
  n <- 150
  G <- matrix(rbinom(n * 5, 3, 0.4), ncol = 5)
  e <- 0.4 * G[, 2] - 0.3 * G[, 4] + rnorm(n)
  got <- fit_variants(G, e)
  ref <- vapply(seq_len(5), function(j) unname(coef(lm(e ~ G[, j]))[2]), numeric(1))
  expect_equal(got$slope, ref, tolerance = 1e-8)
  expect_equal(nrow(got), 5L)
})

test_that("fit_variants returns NA for a monomorphic variant", {
  set.seed(12)
  n <- 100
  got <- fit_variants(cbind(rbinom(n, 3, 0.3), rep(2, n)), rnorm(n))
  expect_false(is.na(got$slope[1]))
  expect_true(is.na(got$slope[2]))
  expect_true(is.na(got$p[2]))
})

test_that("perm_min_p is reproducible and one value per permutation", {
  set.seed(13)
  n <- 120
  G <- matrix(rbinom(n * 8, 3, 0.35), ncol = 8)
  e <- rnorm(n)
  a <- perm_min_p(G, e, n_perm = 50, seed = 99)
  expect_equal(a, perm_min_p(G, e, n_perm = 50, seed = 99))
  expect_length(a, 50L)
  expect_true(all(a >= 0 & a <= 1))
})

test_that("perm_min_p null is not concentrated near zero without signal", {
  set.seed(14)
  n <- 200
  G <- matrix(rbinom(n * 3, 3, 0.4), ncol = 3)
  expect_gt(median(perm_min_p(G, rnorm(n), n_perm = 400, seed = 7)), 0.05)
})

test_that("gene_level_p never returns zero and shrinks with strong signal", {
  expect_equal(gene_level_p(0.5, c(0.1, 0.2, 0.3)), 1 / 4)
  expect_equal(gene_level_p(0.001, rep(0.5, 99)), 1 / 100)
  expect_gt(gene_level_p(0.9, rep(0.5, 99)), 0.9)
})

test_that("a real eQTL beats its own permutation null", {
  set.seed(15)
  n <- 250
  G <- matrix(rbinom(n * 10, 3, 0.3), ncol = 10)
  e <- 0.9 * G[, 3] + rnorm(n)
  obs <- min(fit_variants(G, e)$p, na.rm = TRUE)
  expect_lt(gene_level_p(obs, perm_min_p(G, e, n_perm = 200, seed = 5)), 0.05)
})

test_that("expr_of_gene aligns counts to subject ids and log-transforms", {
  counts <- data.table::data.table(Gene_name = "GENEA", LAB1 = 3, LAB2 = 7)
  meta   <- data.table::data.table(RecordID = c("S1", "S2"), LabID = c("LAB1", "LAB2"))
  expect_equal(expr_of_gene("GENEA", c("S1", "S2"), counts, meta), log2(c(3, 7) + 1))
  expect_true(all(is.na(expr_of_gene("MISSING", "S1", counts, meta))))
})
```

- [ ] **Step 2: Run to verify failure**

Run: `Rscript tests/testthat.R`
Expected: FAIL — `could not find function "fit_variants"`.

- [ ] **Step 3: Write the library**

Create `scripts/lib/eqtl_fit.R`:

```r
# eqtl_fit.R
#
# Vectorized per-variant expression ~ dosage regressions, and the permutation
# machinery for gene-level significance.
#
# Closed-form matrix algebra rather than repeated lm(): the gene-level
# permutation needs n_variants x n_perm fits per gene, far too many for a loop.
# For centred g and e the slope is sum(g*e)/sum(g^2), so all variants and all
# permutations reduce to one matrix product.

center_cols <- function(M) sweep(M, 2, colMeans(M), "-")

#' Fit expression ~ dosage separately for every variant (column) of G.
#' @param G n x m genotype matrix (alt dosage; 0-3 under trisomy).
#' @return data.frame(slope, se, t, p), m rows. Monomorphic variants give NA.
fit_variants <- function(G, e) {
  stopifnot(is.matrix(G), is.numeric(e), nrow(G) == length(e))
  n <- nrow(G)
  if (n < 3) stop("need at least 3 samples")
  Gc  <- center_cols(G)
  ec  <- e - mean(e)
  Sgg <- colSums(Gc^2)
  Sge <- as.vector(crossprod(Gc, ec))
  See <- sum(ec^2)
  slope <- ifelse(Sgg > 0, Sge / Sgg, NA_real_)
  rss   <- See - ifelse(Sgg > 0, Sge^2 / Sgg, 0)
  se    <- ifelse(Sgg > 0, sqrt((rss / (n - 2)) / Sgg), NA_real_)
  tval  <- slope / se
  data.frame(slope = slope, se = se, t = tval, p = 2 * pt(-abs(tval), df = n - 2))
}

#' Minimum p across variants under permutations of expression.
#'
#' Permuting expression breaks the genotype-expression link while preserving
#' genotype LD and the number of variants tested - exactly the multiplicity the
#' "any variant supports the gene" rule ignores.
perm_min_p <- function(G, e, n_perm = 1000, seed = 42) {
  stopifnot(is.matrix(G), nrow(G) == length(e))
  set.seed(seed)
  n <- nrow(G)
  Gc   <- center_cols(G)
  Sgg  <- colSums(Gc^2)
  keep <- Sgg > 0
  if (!any(keep)) return(rep(NA_real_, n_perm))
  Gc  <- Gc[, keep, drop = FALSE]
  Sgg <- Sgg[keep]
  # One n x n_perm matrix of permuted, centred expression; a single matrix
  # product then gives every variant x permutation slope at once.
  E   <- matrix(e[as.vector(replicate(n_perm, sample.int(n)))], nrow = n)
  Ec  <- sweep(E, 2, colMeans(E), "-")
  Sge <- crossprod(Gc, Ec)
  See <- matrix(colSums(Ec^2), nrow = nrow(Sge), ncol = n_perm, byrow = TRUE)
  slope <- Sge / Sgg
  se    <- sqrt(((See - Sge^2 / Sgg) / (n - 2)) / Sgg)
  apply(2 * pt(-abs(slope / se), df = n - 2), 2, min, na.rm = TRUE)
}

#' Gene-level permutation p-value. The +1s keep it strictly positive, which
#' BH-FDR requires.
gene_level_p <- function(min_p_obs, min_p_perm) {
  mp <- min_p_perm[!is.na(min_p_perm)]
  (1 + sum(mp <= min_p_obs)) / (length(mp) + 1)
}

#' Expression vector for one gene, aligned to subject IDs. Counts are keyed by
#' LabID while genotypes are keyed by subject_id, so this maps through metadata.
expr_of_gene <- function(gene_name, subject_ids, counts, meta_t21) {
  row <- counts[Gene_name == gene_name]
  if (nrow(row) == 0) return(rep(NA_real_, length(subject_ids)))
  lab_for_subj <- setNames(meta_t21$LabID, meta_t21$RecordID)
  labs <- lab_for_subj[as.character(subject_ids)]
  log2(suppressWarnings(as.numeric(row[1, match(labs, names(row)), with = FALSE])) + 1)
}
```

- [ ] **Step 4: Run to verify pass**

Run: `Rscript tests/testthat.R`
Expected: PASS.

- [ ] **Step 5: Refactor script 03's regression block**

Add `source("scripts/lib/eqtl_fit.R")` near the top. Inside the `fit_table <- geno_expr[, {...}]` block, replace the `lm()`/`coef()` extraction with:

```r
  fit <- fit_variants(matrix(alt_dosage, ncol = 1), expression)
  .(t21_n     = .N,
    t21_slope = fit$slope[1],
    t21_se    = fit$se[1],
    t21_t     = fit$t[1],
    t21_p     = fit$p[1])
```

- [ ] **Step 6: Verify the refactor is behaviour-preserving**

```bash
cp results/tables/t21_dosage_per_variant.csv /tmp/t21_before.csv
Rscript scripts/03_t21_dosage_boxplots.R
Rscript -e '
suppressPackageStartupMessages(library(data.table))
a <- fread("/tmp/t21_before.csv"); b <- fread("results/tables/t21_dosage_per_variant.csv")
setkey(a, variant_id, ensembl_stable); setkey(b, variant_id, ensembl_stable)
stopifnot(nrow(a) == nrow(b))
cat("max |slope diff|:", max(abs(a$t21_slope - b$t21_slope), na.rm = TRUE), "\n")
cat("max |p diff|    :", max(abs(a$t21_p - b$t21_p), na.rm = TRUE), "\n")
'
```

Expected: both below 1e-8. If not, the refactor changed behaviour — fix before continuing.

- [ ] **Step 7: Append the CHANGELOG block to script 03**

```r
# =============================================================================
# CHANGELOG
# =============================================================================
# 2026-08-31  REPLACED the per-(variant, gene) lm() calls with the vectorized
#             closed-form fits in scripts/lib/eqtl_fit.R. Verified numerically
#             identical to the previous output (max |slope diff| < 1e-8).
#             Reason: the gene-level permutation added below needs
#             n_variants x n_perm fits per gene, beyond lm()'s reach at that
#             scale.
#             Spec: docs/METHODS_SPEC_threshold_and_eqtl_controls.md
```

- [ ] **Step 8: Commit**

```bash
git add scripts/lib/eqtl_fit.R tests/testthat/test-eqtl-fit.R scripts/03_t21_dosage_boxplots.R
git commit -m "Add vectorized eQTL regression library and refactor script 03 onto it"
```

---

## Task 6: Negative controls C1 and C2 inside script 03

**Files:**
- Modify: `scripts/03_t21_dosage_boxplots.R` (new section after the fit table)

**Interfaces:**
- Consumes: `fit_variants()`, `expr_of_gene()` from Task 5.
- Produces: `results/tables/eqtl_negative_controls.csv` with `control, n_genes_tested, n_explained, pct_explained`; rows `observed`, `direction_flip`, `genotype_permutation`.

- [ ] **Step 1: Add the controls section to script 03**

After the fit table is written:

```r
# =============================================================================
# NEGATIVE CONTROLS for the "explained" call
# =============================================================================
# C2 direction flip       - negate the deviation direction and recount. The rule
#                           only asks that some cis variant point the same way
#                           as the deviation, so flipping measures its
#                           discriminating power.
# C1 genotype permutation - shuffle subject labels and refit. Destroys the
#                           genotype-expression link while preserving the number
#                           of variants per gene and their LD, which is the
#                           multiplicity the "any variant" rule ignores.
# A sound criterion collapses toward 0% under both.

cat("\n=== Negative controls ===\n")
N_PERM_SETS <- as.integer(Sys.getenv("T21_PERM_SETS", "20"))

score_explained <- function(dt, deviation_sign, t21_slope, t21_p) {
  supportive <- sign(dt$gtex_slope) == deviation_sign &
                sign(t21_slope) == sign(dt$gtex_slope) &
                !is.na(t21_p) & t21_p < ALPHA_REPRO
  tapply(supportive, dt$Gene_name, function(x) any(x, na.rm = TRUE))
}

obs  <- score_explained(fit_table, fit_table$observed_direction,
                        fit_table$t21_slope, fit_table$t21_p)
flip <- score_explained(fit_table, -fit_table$observed_direction,
                        fit_table$t21_slope, fit_table$t21_p)
cat(sprintf("  observed:          %d/%d explained (%.1f%%)\n",
            sum(obs), length(obs), 100 * mean(obs)))
cat(sprintf("  direction flipped: %d/%d explained (%.1f%%)\n",
            sum(flip), length(flip), 100 * mean(flip)))

gwide <- dcast(geno_t21, subject_id ~ variant_id, value.var = "alt_dosage")
subj  <- gwide$subject_id
G_all <- as.matrix(gwide[, -1])
genes <- unique(fit_table$Gene_name)

perm_rate <- vapply(seq_len(N_PERM_SETS), function(b) {
  set.seed(1000 + b)
  ord <- sample.int(length(subj))
  mean(vapply(genes, function(g) {
    vg   <- fit_table[Gene_name == g]
    cols <- intersect(vg$variant_id, colnames(G_all))
    if (length(cols) == 0) return(FALSE)
    e <- expr_of_gene(g, subj, counts, meta_t21)
    if (all(is.na(e))) return(FALSE)
    fits <- fit_variants(G_all[ord, cols, drop = FALSE], e)
    vgm  <- vg[match(cols, vg$variant_id)]
    any(sign(vgm$gtex_slope) == vgm$observed_direction &
        sign(fits$slope) == sign(vgm$gtex_slope) &
        !is.na(fits$p) & fits$p < ALPHA_REPRO, na.rm = TRUE)
  }, logical(1)))
}, numeric(1))
cat(sprintf("  genotype permuted: %.1f%% explained (mean of %d shuffles, range %.1f-%.1f%%)\n",
            100 * mean(perm_rate), N_PERM_SETS,
            100 * min(perm_rate), 100 * max(perm_rate)))

neg_ctrl <- data.table(
  control        = c("observed", "direction_flip", "genotype_permutation"),
  n_genes_tested = c(length(obs), length(flip), length(genes)),
  n_explained    = c(sum(obs), sum(flip), round(mean(perm_rate) * length(genes))),
  pct_explained  = c(100 * mean(obs), 100 * mean(flip), 100 * mean(perm_rate)))
fwrite(neg_ctrl, "results/tables/eqtl_negative_controls.csv")
if (!file.exists("results/tables/eqtl_negative_controls.csv")) {
  stop("failed to write results/tables/eqtl_negative_controls.csv")
}
cat("  Wrote results/tables/eqtl_negative_controls.csv\n")
```

- [ ] **Step 2: Run and record the rates**

Run: `Rscript scripts/03_t21_dosage_boxplots.R 2>&1 | tail -20`
Expected: three rates. High permuted/flipped rates are the finding, not a script failure.

- [ ] **Step 3: Extend script 03's CHANGELOG**

```r
# 2026-08-31  ADDED negative controls C1 (genotype permutation) and C2
#             (direction flip) -> results/tables/eqtl_negative_controls.csv.
#             Reason: the "explained" call had no null. GTEx and within-T21
#             slopes agree in direction for 93.9% of variants (99.7% at
#             p < 0.05), so the reproducibility criterion is nearly free and the
#             call reduces to a 1-in-N direction coin flip.
```

- [ ] **Step 4: Commit**

```bash
git add scripts/03_t21_dosage_boxplots.R
git commit -m "Add genotype-permutation and direction-flip negative controls to script 03"
```

---

## Task 7: Gene-level permutation p-values in script 03, consumed by script 04

**Files:**
- Modify: `scripts/03_t21_dosage_boxplots.R` (new section), `scripts/04_chr21_lane_assignment.R` (`eqtl_lane` definition)

**Interfaces:**
- Consumes: `fit_variants()`, `perm_min_p()`, `gene_level_p()`, `expr_of_gene()` from Task 5; `genes`, `G_all`, `subj` from Task 6's section.
- Produces: `results/tables/eqtl_gene_level_perm.csv` with `Gene_name, n_variants, min_p_obs, best_variant, p_gene_perm, q_gene_bh, explained_perm`.

- [ ] **Step 1: Add the gene-level permutation section to script 03**

```r
# =============================================================================
# GENE-LEVEL PERMUTATION SIGNIFICANCE (GTEx / FastQTL eGene procedure)
# =============================================================================
# The "any of N variants" rule mostly measures N (21 to 1083 here). Permuting
# expression labels and taking the smallest p across variants per permutation
# builds the null of the BEST variant, handling multiplicity and LD together.
# Unlike picking the lead variant this does not assume the top association is
# causal - in LD the lead variant is frequently only a tag.

N_PERM   <- as.integer(Sys.getenv("T21_N_PERM", "100"))   # 1000 on Alpine
FDR_GENE <- 0.05
cat(sprintf("\n=== Gene-level permutation (%d permutations) ===\n", N_PERM))

perm_res <- rbindlist(lapply(genes, function(g) {
  vg   <- fit_table[Gene_name == g]
  cols <- intersect(vg$variant_id, colnames(G_all))
  e    <- if (length(cols)) expr_of_gene(g, subj, counts, meta_t21) else NA_real_
  if (length(cols) == 0 || all(is.na(e))) {
    return(data.table(Gene_name = g, n_variants = length(cols), min_p_obs = NA_real_,
                      best_variant = NA_character_, p_gene_perm = NA_real_))
  }
  G    <- G_all[, cols, drop = FALSE]
  fits <- fit_variants(G, e)
  best <- which.min(fits$p)
  data.table(Gene_name = g, n_variants = length(cols),
             min_p_obs = fits$p[best], best_variant = cols[best],
             p_gene_perm = gene_level_p(fits$p[best],
                                        perm_min_p(G, e, n_perm = N_PERM, seed = 2026)))
}))
perm_res[, q_gene_bh := p.adjust(p_gene_perm, "BH")]
perm_res[, explained_perm := !is.na(q_gene_bh) & q_gene_bh < FDR_GENE]
setorder(perm_res, p_gene_perm)
print(as.data.frame(perm_res))

fwrite(perm_res, "results/tables/eqtl_gene_level_perm.csv")
if (!file.exists("results/tables/eqtl_gene_level_perm.csv")) {
  stop("failed to write results/tables/eqtl_gene_level_perm.csv")
}
cat("  Wrote results/tables/eqtl_gene_level_perm.csv\n")
```

- [ ] **Step 2: Redefine `eqtl_lane` in script 04**

Replace the `eqtl_lane` `fcase` block with:

```r
# eQTL lane, now gated on gene-level permutation significance rather than
# "at least one supportive variant". The old rule scaled with the number of cis
# variants tested (median n_cis 107 for explained genes vs 36 for the one
# unexplained gene), so it measured variant count more than genetic evidence.
perm <- if (file.exists("results/tables/eqtl_gene_level_perm.csv")) {
  fread("results/tables/eqtl_gene_level_perm.csv")
} else {
  stop("run scripts/03_t21_dosage_boxplots.R first - eqtl_gene_level_perm.csv is missing")
}
m <- merge(m, perm[, .(Gene_name, p_gene_perm, q_gene_bh, explained_perm, best_variant)],
           by = "Gene_name", all.x = TRUE)

m[, eqtl_lane := fcase(
  !(sig_lane %in% c("DE_low", "DE_high")),         "not_evaluated",
  is.na(n_cis_total) | n_cis_total == 0,           "no_GTEx_data",
  !is.na(explained_perm) & explained_perm == TRUE, "explained",
  default =                                        "unexplained")]
```

- [ ] **Step 3: Run 03 then 04 and compare rules**

```bash
Rscript scripts/03_t21_dosage_boxplots.R && Rscript scripts/04_chr21_lane_assignment.R
Rscript -e '
suppressPackageStartupMessages(library(data.table))
p <- fread("results/tables/eqtl_gene_level_perm.csv")
print(p[order(-n_variants), .(Gene_name, n_variants, min_p_obs, p_gene_perm,
                              q_gene_bh, explained_perm)])
'
```

Expected: genes whose old call rested on many variants lose it. Record which.

- [ ] **Step 4: Extend the CHANGELOG blocks**

Script 03:

```r
# 2026-08-31  ADDED gene-level permutation p-values (GTEx/FastQTL eGene
#             procedure) -> results/tables/eqtl_gene_level_perm.csv.
```

Script 04:

```r
# 2026-08-31  REPLACED the eqtl_lane rule "at least one cis variant is
#             direction-matched and reproduces at t21_p < 0.05" with gene-level
#             permutation significance at BH FDR < 0.05.
#             Reason: the old rule scaled with the number of cis variants tested
#             (21 to 1083 per gene, no multiplicity control); median n_cis was
#             107 for explained genes vs 36 for the one unexplained gene, and
#             RBM11 was called explained on 1 supporting variant of 83 where
#             chance predicts ~4. The lead variant was rejected as an
#             alternative because in LD it is frequently a tag, not the causal
#             variant.
```

- [ ] **Step 5: Commit**

```bash
git add scripts/03_t21_dosage_boxplots.R scripts/04_chr21_lane_assignment.R
git commit -m "Gate eqtl_lane on GTEx-style gene-level permutation significance"
```

---

## Task 8: Control gene sets C3 and C5 in scripts 02 and 03

**Files:**
- Modify: `scripts/02_filter_genotypes.R` (control-set selection + genotype stream), `scripts/03_t21_dosage_boxplots.R` (scoring)

**Interfaces:**
- Consumes: `eligible` and `gtex` from Task 4's section of script 02 — **Task 4 must land first**.
- Produces: `data/processed/control_target_variants.csv`, `data/processed/control_genotypes_filtered.csv`, and `results/tables/eqtl_control_geneset_summary.csv` with `gene_set, n_genes, n_with_gtex, n_explained_anyvariant, n_explained_perm, pct_explained_anyvariant, pct_explained_perm`.

- [ ] **Step 1: Add the control gene-set selection to script 02**

After the existing target-variant pull:

```r
# =============================================================================
# CONTROL GENE SETS for the eQTL criterion
# =============================================================================
# C3 Expected_dosage (negative). These sit at the ploidy expectation, so their
#    deviation direction is noise and a direction match is a coin flip. A high
#    "explained" rate here means the criterion is uninformative.
# C5 Strong chr21 eGenes outside the DE set (positive). Genes with the most
#    significant GTEx whole-blood cis-eQTLs should be called explained; if not,
#    the test lacks sensitivity rather than the genes lacking eQTLs.

N_CONTROL <- as.integer(Sys.getenv("T21_N_CONTROL", "30"))
set.seed(2026)

c3 <- eligible[!(Gene_name %in% target_genes$Gene_name) &
                 !is.na(q_outlier) & q_outlier >= OUTLIER_FDR]
c3 <- c3[sample(.N, min(N_CONTROL, .N))]
cat(sprintf("  C3 Expected_dosage control: %d genes\n", nrow(c3)))

best_per_gene <- gtex[!is.na(pval_nominal), .(best_p = min(pval_nominal)),
                      by = ensembl_stable]
c5_stable <- best_per_gene[order(best_p)][
  !ensembl_stable %in% target_genes$ensembl_stable][
  seq_len(min(N_CONTROL, .N)), ensembl_stable]
c5 <- chr21[sub("\\..*$", "", EnsemblID) %in% c5_stable]
cat(sprintf("  C5 strong-eGene control: %d genes\n", nrow(c5)))

control_genes <- rbindlist(list(
  cbind(c3[, .(EnsemblID, Gene_name, norm_log2FC)], gene_set = "expected_dosage"),
  cbind(c5[, .(EnsemblID, Gene_name, norm_log2FC)], gene_set = "strong_egene")),
  fill = TRUE)
control_genes[, ensembl_stable := sub("\\..*$", "", EnsemblID)]
control_genes[, observed_direction := sign(norm_log2FC)]

control_variants <- gtex[ensembl_stable %in% control_genes$ensembl_stable &
                           startsWith(variant_id, "chr21_") &
                           !is.na(pval_nominal) & pval_nominal <= GTEX_PVAL_KEEP]
parsed_c <- tstrsplit(control_variants$variant_id, "_", fixed = TRUE)
control_variants[, POS := as.integer(parsed_c[[2]])]
control_variants <- merge(control_variants,
                          control_genes[, .(ensembl_stable, Gene_name, gene_set,
                                            observed_direction)],
                          by = "ensembl_stable", allow.cartesian = TRUE)
fwrite(control_variants, "data/processed/control_target_variants.csv")
cat(sprintf("  control cis variants: %d (unique positions %d)\n",
            nrow(control_variants), uniqueN(control_variants$POS)))
```

- [ ] **Step 2: Stream genotypes for the control positions**

Reuse script 02's **existing** awk-streaming block, pointed at the control positions and writing `data/processed/control_genotypes_filtered.csv`. Follow the pattern already in the script for the target positions — do not introduce a second mechanism.

- [ ] **Step 3: Run script 02 and confirm both sets**

Run: `Rscript scripts/02_filter_genotypes.R 2>&1 | tail -25`
Expected: both control sets sized, variants written, genotype stream completes. Allow ~5 minutes for the extra 6 GB pass.

- [ ] **Step 4: Score the control sets in script 03**

After the gene-level permutation section, apply **both** rules to the control sets — the any-variant rule from Task 6's `score_explained` and the permutation rule from Task 7 — reading `control_target_variants.csv` and `control_genotypes_filtered.csv`, and write `results/tables/eqtl_control_geneset_summary.csv`. Reporting both side by side is the point: it quantifies how much the permutation rule tightens the false-positive rate.

- [ ] **Step 5: Extend script 02's CHANGELOG**

```r
# 2026-08-31  ADDED control gene sets C3 (Expected_dosage, negative) and C5
#             (strong chr21 eGenes outside the DE set, positive), with their cis
#             variants and genotypes, so the eQTL criterion has a measurable
#             false-positive and sensitivity rate.
```

- [ ] **Step 6: Commit**

```bash
git add scripts/02_filter_genotypes.R scripts/03_t21_dosage_boxplots.R
git commit -m "Add Expected_dosage and strong-eGene control sets for the eQTL criterion"
```

---

## Task 9: Quantitative eQTL prediction in script 04

**Files:**
- Modify: `scripts/04_chr21_lane_assignment.R` (new section before the summary write)

**Interfaces:**
- Consumes: `perm` (Task 7), `data/processed/genotypes_filtered.csv`, `results/tables/t21_dosage_per_variant.csv`, the lane table `m`.
- Produces: `results/tables/eqtl_quantitative_prediction.csv` and `results/figures/eqtl_predicted_vs_observed.{pdf,png}`.

- [ ] **Step 1: Add the quantitative section to script 04**

```r
# =============================================================================
# QUANTITATIVE eQTL PREDICTION
# =============================================================================
# Sign agreement is nearly free (GTEx and within-T21 slopes agree for 93.9% of
# variants), so a directional criterion cannot separate a real explanation from
# a coincidence. This asks whether the eQTL effect size, applied to how far the
# gene's mean alt dosage departs from a diploid-equivalent expectation,
# reproduces the observed deviation:
#
#   predicted_deviation = t21_slope * (mean_dosage - 2 * af),  af = mean_dosage/3
#
# Slope near 1 with high R^2 across genes means the eQTL model explains the
# deviations; slope near 0 means it does not, whatever the sign agreement says.

geno_q <- fread("data/processed/genotypes_filtered.csv")[karyotype == "T21"]
dose   <- geno_q[, .(mean_dosage = mean(alt_dosage, na.rm = TRUE)), by = variant_id]
dose[, af := mean_dosage / 3]
dose[, dosage_excess := mean_dosage - 2 * af]

per_var_q <- fread("results/tables/t21_dosage_per_variant.csv")
quant <- merge(perm[!is.na(best_variant), .(Gene_name, best_variant)],
               per_var_q[, .(Gene_name, variant_id, t21_slope)],
               by.x = c("Gene_name", "best_variant"),
               by.y = c("Gene_name", "variant_id"), all.x = TRUE)
quant <- merge(quant, dose, by.x = "best_variant", by.y = "variant_id", all.x = TRUE)
quant <- merge(quant, m[, .(Gene_name, observed_deviation = norm_log2FC, sig_lane)],
               by = "Gene_name", all.x = TRUE)
quant[, predicted_deviation := t21_slope * dosage_excess]
quant[, residual := observed_deviation - predicted_deviation]

qfit <- lm(observed_deviation ~ predicted_deviation, data = quant)
cat("\n=== Quantitative eQTL prediction ===\n")
print(summary(qfit)$coefficients)
cat(sprintf("  R^2 = %.3f, slope = %.2f\n", summary(qfit)$r.squared, coef(qfit)[2]))

fwrite(quant, "results/tables/eqtl_quantitative_prediction.csv")

qp <- ggplot(quant, aes(x = predicted_deviation, y = observed_deviation)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_smooth(method = "lm", se = TRUE, colour = "#2166AC", linewidth = 0.5) +
  geom_point(aes(colour = sig_lane), size = 2) +
  ggrepel::geom_text_repel(aes(label = Gene_name), size = 2.6,
                           max.overlaps = Inf, seed = 1) +
  labs(x = "Predicted deviation from cis-eQTL (log2)",
       y = "Observed ploidy-corrected deviation (log2)",
       title = "Does the eQTL predict the magnitude of the deviation?",
       subtitle = sprintf("dashed = y = x; R^2 = %.3f, slope = %.2f",
                          summary(qfit)$r.squared, coef(qfit)[2]),
       colour = NULL) +
  theme_bw(base_size = 9)
ggsave("results/figures/eqtl_predicted_vs_observed.pdf", qp, width = 7, height = 5,
       units = "in", device = "pdf")
ggsave("results/figures/eqtl_predicted_vs_observed.png", qp, width = 7, height = 5,
       units = "in", dpi = 300)
for (f in c("results/tables/eqtl_quantitative_prediction.csv",
            "results/figures/eqtl_predicted_vs_observed.pdf",
            "results/figures/eqtl_predicted_vs_observed.png")) {
  if (!file.exists(f)) stop("failed to write ", f)
  cat("  Wrote ", f, "\n", sep = "")
}
```

- [ ] **Step 2: Run and record**

Run: `Rscript scripts/04_chr21_lane_assignment.R 2>&1 | tail -25`
Expected: a printed regression, R^2, slope, three files. Record slope and R^2 — these decide whether "explained by eQTL" survives as a claim.

- [ ] **Step 3: Extend script 04's CHANGELOG**

```r
# 2026-08-31  ADDED a quantitative test of whether the cis-eQTL predicts the
#             MAGNITUDE of the deviation, not only its sign
#             -> results/tables/eqtl_quantitative_prediction.csv.
#             Reason: GTEx and within-T21 slopes agree in direction for 93.9% of
#             variants, so sign agreement cannot separate explanation from
#             coincidence.
```

- [ ] **Step 4: Commit**

```bash
git add scripts/04_chr21_lane_assignment.R
git commit -m "Add quantitative test of whether eQTL effect size predicts deviation magnitude"
```

---

## Task 10 (conditional): Alpine offload for the full permutation run

**Run only if the gate in Step 1 fails.** Permutation work defaults to `T21_N_PERM = 100` locally; the publishable run wants 1000 plus 60 control genes.

**Files:**
- Create: `scripts/alpine/ralpine`, `scripts/alpine/permutation_controls.sbatch`
- Modify: `scripts/download_vcf_array.sh`

- [ ] **Step 1: Gate — measure the local cost**

```bash
Rscript -e '
source("scripts/lib/eqtl_fit.R")
set.seed(1); n <- 302; G <- matrix(rbinom(n*1083, 3, .3), ncol = 1083); e <- rnorm(n)
cat("1000 perms on the largest gene (1083 variants):",
    round(system.time(perm_min_p(G, e, n_perm = 1000, seed = 1))["elapsed"], 1), "s\n")
'
```

If the full run extrapolates to **under 30 minutes**, skip this task and run locally with `T21_N_PERM=1000`. Record the measurement either way.

- [ ] **Step 2: Install the Alpine access boundary**

No absolute paths go into tracked files. The source checkout and the remote
root are both supplied by the environment, so this works on any machine and for
any account.

```bash
mkdir -p scripts/alpine logs
# FM_PDO_REPO: your fm-pdo-evaluator checkout, wherever it lives.
: "${FM_PDO_REPO:?set FM_PDO_REPO to your fm-pdo-evaluator checkout}"
cp "$FM_PDO_REPO/scripts/alpine/ralpine" scripts/alpine/ralpine
chmod +x scripts/alpine/ralpine
```

Then edit `scripts/alpine/ralpine` so the remote root carries **no hardcoded
path**. Replace its `REMOTE_ROOT` assignment with:

```bash
# Remote checkout location. Deliberately has no default: hardcoding a
# /projects/<account>/... path bakes one user's account into version control and
# breaks for everyone else. Export ALPINE_ROOT in your shell profile, e.g.
#   export ALPINE_ROOT="$ALPINE_PROJECTS/repositories/T21-eQTL"
REMOTE_ROOT="${ALPINE_ROOT:?set ALPINE_ROOT to the Alpine checkout path}"
```

Update the header comment to name this repo. Keep the read-only allowlist and
the `submit`/`cancel` restrictions **exactly as they are** — they are a
deliberate safety boundary, not boilerplate.

- [ ] **Step 3: Write the sbatch script**

Create `scripts/alpine/permutation_controls.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=t21-perm
#SBATCH --partition=acpu          # amilan was retired by CURC on 2026-08-24
#SBATCH --qos=cpu-normal          # acpu accepts cpu-normal / cpu-long, not normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# Full-scale permutation run for the eQTL controls: scripts 03 and 04 at
# N_PERM = 1000 over the DE genes and both control gene sets. Script 02 must
# have run already so the genotype and variant tables exist on this checkout.
#
# Submit: scripts/alpine/ralpine submit scripts/alpine/permutation_controls.sbatch
set -euo pipefail

REPO="${REPO:-$SLURM_SUBMIT_DIR}"
cd "$REPO"

# R provisioning on Alpine is UNVERIFIED for this repo. `module load` is known to
# fail inside batch jobs for the anaconda module in fm-pdo-evaluator; whether the
# R module behaves the same here is untested. Try the module, then an optional
# caller-supplied prefix, then whatever is already on PATH - and fail loudly
# rather than running a wrong R.
#
# No absolute paths here: R_PREFIX is supplied by the environment (export it in
# your shell profile or pass it with --export to sbatch). Hardcoding a
# /projects/<account>/... path would bake one user's account into the repo.
if module load "${R_MODULE:-R}" 2>/dev/null && command -v Rscript >/dev/null; then
  echo "R via module ${R_MODULE:-R}"
elif [[ -n "${R_PREFIX:-}" && -x "${R_PREFIX}/bin/Rscript" ]]; then
  export PATH="${R_PREFIX}/bin:$PATH"
  echo "R via R_PREFIX"
elif command -v Rscript >/dev/null; then
  echo "R already on PATH"
else
  echo "ERROR: no usable R. Load an R module or export R_PREFIX to an R install." >&2
  exit 1
fi

Rscript -e 'cat(R.version.string, "\n")'

export T21_N_PERM="${T21_N_PERM:-1000}"
export T21_PERM_SETS="${T21_PERM_SETS:-100}"
export T21_N_CONTROL="${T21_N_CONTROL:-30}"
echo "N_PERM=$T21_N_PERM PERM_SETS=$T21_PERM_SETS N_CONTROL=$T21_N_CONTROL"

Rscript scripts/03_t21_dosage_boxplots.R
Rscript scripts/04_chr21_lane_assignment.R

echo "done: $(date)"
```

- [ ] **Step 4: Fix the retired partition in the VCF download script**

In `scripts/download_vcf_array.sh`, replace:

```bash
#SBATCH --partition=amilan
#SBATCH --qos=normal
```

with:

```bash
#SBATCH --partition=acpu       # amilan was retired by CURC on 2026-08-24
#SBATCH --qos=cpu-normal       # acpu accepts cpu-normal / cpu-long, not normal
```

and append a CHANGELOG block recording the change and its reason.

- [ ] **Step 5: Submit and collect**

Both `ALPINE_ROOT` and, if the R module does not work, `R_PREFIX` must be set in
the environment first — neither has a hardcoded default.

```bash
git push
scripts/alpine/ralpine update
scripts/alpine/ralpine submit scripts/alpine/permutation_controls.sbatch
scripts/alpine/ralpine sq
scripts/alpine/ralpine log t21-perm
scripts/alpine/ralpine pull results/tables/eqtl_gene_level_perm.csv results/tables/
scripts/alpine/ralpine pull results/tables/eqtl_control_geneset_summary.csv results/tables/
```

`ralpine submit` refuses a stale checkout, so `update` must succeed first.

- [ ] **Step 6: Commit**

```bash
git add scripts/alpine/ralpine scripts/alpine/permutation_controls.sbatch scripts/download_vcf_array.sh
git commit -m "Add Alpine offload for the full permutation run; fix retired amilan partition"
```

---

## Task 11: Figure, decision log, and pipeline documentation

**Files:**
- Modify: `scripts/07_three_panel_figure.R`, `docs/REPO_STATE.md`, `CLAUDE.md`
- Delete: `scripts/05_chr21_distribution_panel.R`

- [ ] **Step 1: Regenerate the figure**

Run: `Rscript scripts/07_three_panel_figure.R 2>&1 | tail -10`
Expected: runs clean; labelled gene groups come from the lane table, so the new DE set is picked up automatically.

- [ ] **Step 2: Visually check it**

Open `results/figures/three_panel_summary.png`. Confirm: no overlapping labels, exactly one legend, panel C subtitle reports the new outside-noise count.

- [ ] **Step 3: Append the decision log to `docs/REPO_STATE.md`**

Append at the very bottom:

```markdown
---

## Decision log: 2026-08-31 Hunter compliance, threshold, and eQTL controls

Spec: `docs/METHODS_SPEC_threshold_and_eqtl_controls.md`.
Plan: `docs/superpowers/plans/2026-08-31-threshold-and-eqtl-controls.md`.

**0. The ploidy arm had no library-size normalization.** `01:210` set
`normalizationFactors(dds_norm) <- norm_matrix`, where norm_matrix holds only
1.0 and 1.5. normalizationFactors replace size factors outright and DESeq()
skips estimateSizeFactors() when they are set, so every sample was treated as
equal depth across a 21.1M-52.3M read range (2.47x, CV 12.2%). Replaced with
`estimateSizeFactors(dds_norm, normMatrix = norm_matrix, controlGenes = ...)`,
which is what Hunter et al. specify and what the raw arm already did. This
changed every ploidy-corrected number in the pipeline.

**1. Hunter et al.'s five steps, audited.** Steps 1 and 2 were already
satisfied, implemented as LABELS rather than filters - which is the deliberate
design here, so that genes below threshold stay in the table and remain
available for eQTL-support analysis. Step 1: the absolute count form
(total >= 30) excludes nothing (min 437, median 195,547), and the q20 baseMean
rule (25.10, flags 32 of 160) is a strict subset of Hunter's baseMean < 30
(flags 34) - a 2-gene difference. Step 2: ENCODE hg38-blacklist.v2 region
overlap via process_blacklist.R (10 genes) plus 9 KNOWN_REPEAT_GENES, at gene
level; read-level masking is unavailable because counts arrive precomputed from
Synapse. Step 3 was correct in the raw arm and broken in the ploidy arm; see
above. Step 4 is N/A at n = 399, but its residue (DESeq2 nulling p-values for
count outliers, e.g. MX1) now has a cooksCutoff=FALSE sensitivity arm. Step 5
was correct in intent.

An earlier draft of this plan proposed adding a hard rowSums >= 30 pre-filter.
That was WRONG and was removed: it would have emptied the Low_expression and
High_repeats lanes the Sankey depends on, and discarded genes needed for later
eQTL work.

**2. The cohort-noise null was the wrong reference.** Ploidy normalization acts
only on chr21 (mean |raw - norm| 0.583 on chr21, 0.0048 off it), so a diploid
gene's "corrected" log2FC is its ordinary T21-vs-Control log2FC. Replaced with a
chr21-internal null.

**3. MAD, not SD.** chr21 SD/MAD = 1.49; removing 5 of 160 genes moves SD 18.9%
and MAD 1.1%. Those 5 include OLIG2, so an SD rule let a primary candidate help
set the bar it had to clear.

**4. Estimate the scale after the expression filter.** The naive top hits are
KRTAP12-4 (baseMean 1.4), CLIC6 (5.6), KRTAP10-9 (1.3) - count noise.

**5. k is set by FDR, not by hand.** The chr21 bulk is normal once the 10% most
extreme are trimmed (Shapiro p = 0.17 vs 6e-6 untrimmed). The old k = 1 flagged
36.2% of genes against a 31.7% null expectation - a 1.1x enrichment. Re-derive
k after the size-factor fix.

**6. 90.5% of the chr21 spread was real biology**, not measurement error (RMS
sampling SE 0.108 vs observed SD 0.352) - measured pre-fix; re-measure. The
filter is an outlier-definition question, not a confidence question. DE_high was
stable across thresholds while DE_low ranged 4 to 18 genes; always report the
sensitivity table.

**7. The eQTL "explained" call had no null.** GTEx and within-T21 slopes agree
in direction for 93.9% of variants (99.7% at p<0.05), so the reproducibility
criterion is nearly free and the call reduces to a 1-in-N direction coin flip.
Added negative controls (genotype permutation, direction flip, Expected_dosage
set) and positive controls (GTEx-vs-T21 concordance, already passing at Spearman
0.841; strong eGenes).

**8. Gene-level permutation, not the lead variant.** In LD the top-associated
variant is frequently a tag rather than the causal one, so substituting it
trades one arbitrary rule for another. The GTEx/FastQTL eGene permutation
handles multiplicity and LD together.

**9. Quantitative prediction is required.** Sign agreement cannot separate
explanation from coincidence; the eQTL must predict the magnitude.

**10. The analysis cohort is asymmetric: 302 T21 + 95 Control = 397.** T21
require RNA-seq and chr21 WGS (302 of 304); Controls require RNA-seq only (95 of
95), because genotypes are used solely for the within-T21 dosage regressions
that controls never enter - requiring WGS of controls would discard 89 of 95 for
no analytic gain. DE now runs on the same 397 people as the eQTL step, removing
the old "302 of 304" mismatch. Only 6 controls have genotypes at all, so there
is no usable internal diploid arm for eQTL work; GTEx is the only diploid
reference. No race or ethnicity is recorded in the available metadata.

**Deliberately NOT done in this pass:** headline numbers and figures in
`CLAUDE.md` and `README.md` were left untouched pending settled results. Only
the pipeline description was updated.
```

- [ ] **Step 4: Update the pipeline description in CLAUDE.md — no numbers**

Update **only**:
- the pipeline table: script 00 now defines the analysis cohort and emits Table 1; scripts 01-04 carry the cohort restriction, the size-factor fix, the FDR threshold, the eQTL controls, the gene-level permutation, and the quantitative test, with their new output files;
- the methodology section: replace the "Cohort-noise (within-cohort SD) filter" description with the chr21-internal median/MAD + FDR rule, replace the eQTL "supported" definition with the permutation-based one, and record the Hunter et al. compliance audit including the step 2 limitation;
- the constants list: replace `MAGNITUDE_THRESHOLD = 1.0` with `OUTLIER_FDR = 0.10`, add `N_PERM`, `FDR_GENE` (and `LOW_EXPR_ABS` if the low-expression label was switched to the absolute form); state the label-not-filter principle;
- the input list: correct script 01's chr21 output filename to `deseq2_chr21_genes_both_analyses.csv`;
- the Quick Start block: remove the deleted `05_chr21_distribution_panel.R`;
- the repository layout: add `scripts/lib/` and `scripts/alpine/`.

**Do NOT touch** the "Headline result" block, any reported counts, or any figure description. Leave `README.md` entirely alone.

- [ ] **Step 5: Delete the superseded script**

```bash
git rm scripts/05_chr21_distribution_panel.R
```

- [ ] **Step 6: Run the full pipeline from a clean state**

```bash
Rscript scripts/01_deseq2_analysis.R && \
Rscript scripts/02_filter_genotypes.R && \
Rscript scripts/03_t21_dosage_boxplots.R && \
Rscript scripts/04_chr21_lane_assignment.R && \
Rscript scripts/05_alluvial_lane_assignment.R && \
Rscript scripts/06_chr21_distribution_panel.R && \
Rscript scripts/07_three_panel_figure.R
```

Expected: every script exits 0. This is the check that no script reads a file no script writes — the failure mode behind the original stale-input bug.

- [ ] **Step 7: Run the unit tests**

Run: `Rscript tests/testthat.R`
Expected: PASS.

- [ ] **Step 8: Commit and open the PR**

```bash
git add CLAUDE.md scripts/07_three_panel_figure.R
git commit -m "Document pipeline changes; regenerate summary figure"
git push -u origin hunter-threshold-eqtl-controls
gh pr create --title "Hunter compliance, robust chr21 threshold, and eQTL controls" \
  --body "$(cat <<'EOF'
## Summary
- Fixes the ploidy arm's size-factor calculation. It ran with NO library-size
  normalization because normalizationFactors were assigned directly; every
  ploidy-corrected number in the pipeline changes.
- Restricts DE to the analysis cohort: T21 need RNA-seq and WGS (302 of 304),
  Controls need RNA-seq only (95 of 95), so DE and eQTL cover the same 397
  people. Adds Table 1.
- Audits Hunter et al.'s five steps: steps 1 and 2 were already satisfied as
  labels (not filters, by design); adds a Cook's-distance sensitivity arm for
  step 4 and documents the repeat-masking limitation.
- Replaces the cohort-SD magnitude filter with an FDR-controlled robust outlier
  test against a chr21-internal median/MAD null.
- Adds five eQTL controls; replaces the any-variant "explained" rule with
  GTEx-style gene-level permutation p-values plus BH-FDR; adds a quantitative
  test of whether eQTL effect size predicts deviation magnitude.
- Adds Table 1 for the RNA-seq + WGS analysis cohort (302 T21, 6 Control).
- All changes integrated into existing pipeline stages; each modified script
  carries a CHANGELOG block recording what was replaced and why.

Headline numbers in CLAUDE.md and README.md are deliberately NOT updated in this
PR; only the pipeline description is.

## Spec
docs/METHODS_SPEC_threshold_and_eqtl_controls.md

## Test plan
- `Rscript tests/testthat.R` passes.
- Full 01-07 chain runs clean from a clean state.
EOF
)"
```

---

## Self-Review

**Spec coverage:**

| Spec section | Task |
|---|---|
| Hunter step 1 (coverage filter) | 1 |
| Hunter step 2 (repeat masking) | 1 (documented limitation) |
| Hunter step 3 (chr21 out of size factors) | 1 |
| Hunter step 4 (noisy samples / Cook's) | 1 |
| Hunter step 5 (ploidy null) | 1 |
| Problem 1 (wrong null reference) | 2, 3 |
| Problem 2 (MAD not SD) | 2, 3 |
| Problem 3 (filter before scale) | 2 documents, 3 implements |
| Problem 4 (k by FDR; sensitivity table) | 2, 3 |
| Problem 5 (MX1 / Cook's) | 1 |
| Problem 6 C1 genotype permutation | 5 |
| Problem 6 C2 direction flip | 5 |
| Problem 6 C3 Expected_dosage set | 7 |
| Problem 6 C4 GTEx-vs-T21 concordance | passing already; logged in Task 11 Step 3 |
| Problem 6 C5 strong eGenes | 7 |
| Fix: gene-level permutation, not lead variant | 6 |
| Fix: BH-FDR across genes | 6 |
| Fix: quantitative prediction | 8 |
| Problem 7 (Table 1) | 9 |
| Known data limitation (6 controls) | 9 reports; Task 11 Step 3 logs |

No gaps.

**Placeholder scan:** no TBD/TODO; every code step carries runnable code; no "similar to Task N" back-references. Task 8 Step 2 and Task 11 Step 4 describe reusing an existing in-script mechanism and editing prose — both name the exact target and constraint.

**Type consistency:** `chr21_null()` returns `list(center, scale, n)`, consumed as `null$center`/`null$scale` in Tasks 3, 4. `fit_variants()` returns `data.frame(slope, se, t, p)`, consumed as `fit$slope[1]`, `fits$p`, `fits$slope` in Tasks 5-8. `perm_min_p()` returns numeric, consumed by `gene_level_p()` in Tasks 5, 7, 10. `expr_of_gene(gene_name, subject_ids, counts, meta_t21)` has one signature at its definition (Task 5) and all call sites (Tasks 6, 7). Columns `dev_z`/`q_outlier` written in Task 4, read in Tasks 4, 8. `best_variant` produced in Task 7, consumed in Task 9. Table 1 helpers are used only in Task 1.

**Known risks to watch during execution:**
1. **Task 2 invalidates every pre-measured number in this plan.** The chr21 null values, the k table, the DE gene lists, and the 90.5%-biology figure were all measured on the unnormalized arm. Re-measure at Task 3 Step 6 and do not carry the old numbers forward.
2. Task 2 adds no gene filter, so `norm_matrix` construction is untouched and its existing `stopifnot` dimension checks still apply. If any future step is tempted to filter genes before DESeq2, stop: it would empty the `Low_expression` and `High_repeats` lanes.
3. The genotype join key is `subject_id`, derived from `LabID` by stripping the trailing visit suffix - **not** `RecordID`, which yields zero matches (verified). Task 1 Step 8 asserts 302 / 95 as the guard. Task 6's controls use the same mapping via `expr_of_gene()`.
4. Task 8 Step 1 references `eligible` and `gtex`, defined by Task 4 and existing code in script 02. Executing Task 8 before Task 4 will fail on a missing `eligible`.
5. Task 10's R provisioning on Alpine is **unverified**. The sbatch tries `module load`, falls back to PATH, and fails loudly rather than running a wrong R. Do not assert either path works until a job completes.
