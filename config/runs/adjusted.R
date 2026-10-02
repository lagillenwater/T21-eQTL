# Covariate- and composition-adjusted run; mosaic T21 excluded.
list(
  name        = "adjusted",
  covariates  = c("Age_at_visit", "Sex", "BMI", "Sample_source"),
  composition = list(
    source    = "cytof",
    path      = "data/synapse/syn31488783/HTP_CyTOF_CD45posCD66low_FlowSOM_cluster_percentage_Synapse.txt",
    reference = "Classical monocytes and M-MDSCs"),
  exclude     = list(karyotype_subtype = "mosaic_T21"),
  ploidy      = 1.5,
  thresholds  = list(alpha_de = 0.01, deviation_lfc = log2(1.5),
                     deviation_lfc_t2 = log2(4/3), low_expr_basemean = 30,
                     fdr_gene = 0.05, gtex_pval_keep = 1e-4,
                     outlier_fdr = 0.10, alpha_repro = 0.05,
                     # standalone eQTL controls (script 02 selects, 03 tests)
                     n_positive_controls = 10, decoy_min_distance = 5e6,
                     # controls v2: matched one-per-locus positives (GTEx eGene
                     # q, TSS separation); k decoy sets per gene, MAF floor
                     positive_egene_qval = 1e-4, positive_min_separation = 1e6,
                     positive_dev_separation = 1e5, positive_min_variants = 10,
                     n_decoy_sets = 5, decoy_min_maf = 0.05)
)
