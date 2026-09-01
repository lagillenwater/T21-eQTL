# Tight plan: Hunter et al. extended to the HTP cohort

Supersedes Tasks 8-11 of `2026-08-31-threshold-and-eqtl-controls.md`.
Tasks 1-7 of that plan stand (cohort, DESeq2 fix, threshold library,
eQTL library, negative controls, gene-level permutation).

**Question.** Do chr21 genes in T21 whole blood follow the 1.5x dosage
expectation, and where they don't, is there a cis-eQTL?

| step | what | why | status |
|---|---|---|---|
| 1 | Cohort: 302 T21 (RNA-seq + WGS), 95 controls. Table 1. | defines the analysis set | done |
| 2 | Trisomy-aware DESeq2 as Hunter: ploidy-normalised, chr21 excluded from size factors | replicates their method | done, bug fixed |
| 3 | Classify chr21 genes as Hunter did: deviating = padj < 0.01 AND abs(corrected log2FC) >= log2(1.5); label low-expression and repeat genes, do not drop them | their FC >= 1.5 cut on the corrected scale; the effect-size floor is the control against large-n over-calling | Task A |
| 4 | Composition control: flag deviating genes whose co-expressed non-chr21 partners also shift | whole blood is not their LCLs; composition can mimic deviation | Task B |
| 5 | Per deviating gene, gene-level permutation test for a cis-eQTL in T21. Report has / has not. | the any-variant rule failed its control; permutation replaces it | done |
| 6 | Report: their finding replicates (median FC 1.47); N deviate; which are programs vs gene-specific; which gene-specific ones carry a cis-eQTL | | Task C |

**Not in this plan:** the FDR-outlier test (stays as annotation columns
`dev_z`, `q_outlier`; does not drive classification) and any "explained by
eQTL" claim.

**One dial:** `DEVIATION_LFC = log2(1.5)`. Hunter's value. Do not move it
after looking at gene names. `log2(4/3)` is the only alternative worth
considering, and it is a one-line change.

---

## Task A: Hunter's classification rule (step 3)

Files: `scripts/02_filter_genotypes.R`, `scripts/04_chr21_lane_assignment.R`.

- Add `DEVIATION_LFC <- log2(1.5)` to both scripts' constants.
- Script 02 target rule becomes
  `!is.na(norm_padj) & norm_padj < ALPHA_DE & abs(norm_log2FC) >= DEVIATION_LFC`
  over the `eligible` table. Remove the `q_outlier` condition from the target
  rule; leave `dev_z`/`q_outlier` computed as annotations.
- Script 04: `passes_magnitude_filter := abs(norm_log2FC) >= DEVIATION_LFC`
  for eligible genes (NA for ineligible). The `sig_lane` fcase order stays
  exactly as it is (high_repeat, low_expr, then the filter).
- Keep `OUTLIER_FDR`, `dev_z`, `q_outlier`, `chr21_k_sensitivity.csv` as
  annotations. Do not delete them.
- Re-run 02 -> 03 -> 04. Expected deviating set: OLIG2, COL6A1, RIPK4,
  TSPEAR (4). Lane table stays 160 rows.
- CHANGELOG entries in both scripts. Tests still pass.

## Task B: composition control (step 4)

Files: create `scripts/lib/composition.R` + `tests/testthat/test-composition.R`;
modify `scripts/04_chr21_lane_assignment.R`.

- `partner_shift(gene, L_ctrl, L_all_gene_row, lfc, n_partners = 20)`:
  correlate the gene's log2-CPM (controls only) with every expressed non-chr21
  gene; take the top `n_partners`; return their median T21-vs-Control log2FC.
- `partner_null(lfc_bg, n_partners, n_draw = 2000, seed)`: median log2FC of
  `n_draw` random `n_partners`-gene sets. Vectorised, no loops over draws.
- `composition_verdict(gene_lfc, partner_lfc, p)`: PROGRAM if p < 0.05 and
  partner_lfc / gene_lfc >= 0.5; MIXED if p < 0.05; else GENE-SPECIFIC.
- Script 04 runs this for every gene in a deviating lane and writes
  `results/tables/chr21_composition_control.csv`:
  `Gene_name, gene_lfc, partner_lfc, p_partners, program_share, residual_lfc,
  verdict`, where `residual_lfc = gene_lfc - partner_lfc`.
- Unit tests: a synthetic gene whose partners all shift -> PROGRAM; whose
  partners sit at zero -> GENE-SPECIFIC; null is reproducible under seed.
- CHANGELOG entry. Tests pass.

## Task C: report and drift check (step 6)

- `docs/REPO_STATE.md`: append the decision log (the one dial, what was
  dropped and why).
- `CLAUDE.md`: pipeline description only, no result numbers.
- Regenerate `scripts/07_three_panel_figure.R`.
- Write `results/tables/hunter_extension_summary.md`: one table per step 6.
- **Drift check:** diff what was built against the six-step table above.
  Report every place the implementation departs from the table.
