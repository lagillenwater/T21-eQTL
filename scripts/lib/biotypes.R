# biotypes.R
#
# The chr21 target gene set, by GENCODE biotype: the biotypes GTEx whole
# blood tests for cis-eQTLs (protein-coding, lncRNA, and the pseudogene
# classes). Every script that restricts a gene table sources this file, so
# the set is defined once. Widened from protein-coding only on 2026-09-04.
# The low-expression floor (Hunter et al.'s baseMean < 30, scripts 02 and 04)
# is absolute, so widening the set does not move it.

TARGET_BIOTYPES <- c("protein_coding", "lncRNA",
                     "processed_pseudogene",
                     "transcribed_processed_pseudogene",
                     "transcribed_unprocessed_pseudogene",
                     "unprocessed_pseudogene", "pseudogene")

# Short label for figure text.
TARGET_BIOTYPES_LABEL <- "protein-coding, lncRNA, pseudogene"
