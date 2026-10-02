# karyotype_subtype.R
#
# Karyotype subtype from the INCLUDE Data Hub export. The subtype file keys
# on INCLUDE participant ids (pt-...); the mapping file gives the HTP external
# participant id (e.g. HTP0382). Labels are MONDO disease codes as exported:
#   T21                complete trisomy 21
#   DS_T21             Down syndrome, subtype unspecified
#   mosaic_T21         mosaic trisomy 21
#   translocation_T21  translocation Down syndrome
#   D21                euploid
# Non-HTP participants (KFDS ids) are dropped.

#' One row per HTP subject: (subject_id, karyotype_subtype).
karyotype_subtype_table <- function(subtype_path, mapping_path) {
  s <- data.table::fread(subtype_path, header = TRUE)
  data.table::setnames(s, c("row", "participant_id", "karyotype_subtype"))
  m <- data.table::fread(mapping_path, header = TRUE)
  data.table::setnames(m, c("row", "participant_id", "subject_id", "htp_num", "participant2"))
  j <- merge(s[, .(participant_id, karyotype_subtype)],
             m[, .(participant_id, subject_id)], by = "participant_id")
  unique(j[grepl("^HTP[0-9]+$", subject_id), .(subject_id, karyotype_subtype)])
}

#' Attach karyotype_subtype to metadata by subject_id; NA when unmapped.
#' Row order and count of `meta` are preserved.
add_karyotype_subtype <- function(meta, subtypes) {
  stopifnot("subject_id" %in% names(meta))
  m <- data.table::as.data.table(meta)
  lookup <- setNames(subtypes$karyotype_subtype, subtypes$subject_id)
  m[, karyotype_subtype := unname(lookup[as.character(subject_id)])]
  m
}
