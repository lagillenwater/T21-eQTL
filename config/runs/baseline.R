# Hunter-style primary run: no covariates, no exclusions.
list(
  name        = "baseline",
  covariates  = character(0),
  composition = NULL,
  exclude     = list(),
  ploidy      = 1.5,
  thresholds  = list(alpha_de = 0.01, deviation_lfc = log2(1.5),
                     deviation_lfc_t2 = log2(4/3), low_expr_basemean = 30,
                     fdr_gene = 0.05, gtex_pval_keep = 1e-4,
                     outlier_fdr = 0.10, alpha_repro = 0.05)
)
