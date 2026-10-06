# The adjusted run with karyotype x composition interactions: every CyTOF
# cell fraction may relate to expression differently in T21 and D21. The
# fractions are centred, so the karyotype effect is estimated at the
# cohort-mean composition. Mosaic T21 excluded, as in adjusted.
list(
  name        = "adjusted_interaction",
  covariates  = c("Age_at_visit", "Sex", "BMI", "Sample_source"),
  composition = list(
    source    = "cytof",
    path      = "data/synapse/syn31488783/HTP_CyTOF_CD45posCD66low_FlowSOM_cluster_percentage_Synapse.txt",
    reference = "Classical monocytes and M-MDSCs",
    interact_karyotype = TRUE),
  exclude     = list(karyotype_subtype = "mosaic_T21"),
  ploidy      = 1.5,
  thresholds  = list(alpha_de = 0.01, deviation_lfc = log2(1.5),
                     deviation_lfc_t2 = log2(4/3), low_expr_basemean = 30,
                     fdr_gene = 0.05, gtex_pval_keep = 1e-4,
                     outlier_fdr = 0.10, alpha_repro = 0.05,
                     distal_min_distance = 5e6, distal_min_maf = 0.05)
)
