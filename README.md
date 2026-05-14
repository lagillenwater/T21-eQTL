# T21-eQTL: Dosage Compensation in Down Syndrome

Cohort-scale extension of Hunter et al. (2023) "Transcription dosage
compensation does not occur in Down syndrome" (BMC Biology 21:228), applied
to the Human Trisome Project (HTP): 304 T21 + 95 Control whole-blood RNA-seq
samples, 302 of the T21 subjects with paired chr21 genotypes.

The pipeline asks, for each chromosome 21 protein-coding gene whose
ploidy-corrected expression deviates from the 1.5x trisomy expectation by
more than typical genome-wide cohort noise: can the deviation be explained
by a common cis-acting eQTL in GTEx whole blood that reproduces in our T21
cohort, or does it remain a candidate for true regulatory dosage
compensation?

## Headline result

Of the 160 chr21 protein-coding genes tested:

| Lane | Count |
|---|---|
| Expected dosage (within cohort noise) | 119 |
| DE_high, eQTL-supported | 6 |
| DE_low, eQTL-supported | 9 |
| DE_low, eQTL-tested but unsupported | 3 (BACE2, ADARB1, NRIP1) |
| DE_low, no GTEx cis-eQTL coverage | 3 (CLIC6, JAM2, FTCD-class) |
| High repeats / Low expression / Not DE | 20 |

All higher-than-expected genes (DE_high) are eQTL-supported - acting as a
positive control. Among lower-than-expected genes that could be tested
against GTEx, 9 of 12 (75%) are eQTL-supported; the remaining 3 (BACE2,
ADARB1, NRIP1) plus the 3 untestable ones are the residual compensation
candidates.

## Quick start

```bash
Rscript install_packages.R                        # one-time
Rscript scripts/00_preprocess_data.R              # long -> wide counts
Rscript scripts/01_deseq2_analysis.R              # trisomy-aware DESeq2
Rscript scripts/02_filter_genotypes.R             # gene + variant universe
Rscript scripts/03_t21_dosage_boxplots.R          # within-T21 fits
Rscript scripts/04_chr21_lane_assignment.R        # MAIN per-gene table
Rscript scripts/05_alluvial_lane_assignment.R     # alluvial + Sankey
Rscript scripts/06_chr21_distribution_panel.R     # population view
```

The full 00-06 chain runs in ~30 minutes. Legacy and supplementary scripts
(paper-style Panel D, additional eQTL diagnostics, supporting figure
panels) live under `scripts/archive/`.

## Documentation

See [CLAUDE.md](CLAUDE.md) for the canonical project documentation:
methodology, pipeline walkthrough, data inputs, full script catalog, output
file inventory, code style, and common gotchas.

## Citation

Hunter, S., Hendrix, J., Freeman, J., Dowell, R.D., & Allen, M.A. (2023).
Transcription dosage compensation does not occur in Down syndrome.
*BMC Biology* 21:228. https://doi.org/10.1186/s12915-023-01700-4

## License

BSD 3-Clause License — see [LICENSE](LICENSE).
