test_that("karyotype_subtype_table joins INCLUDE ids to HTP subjects and keeps HTP rows only", {
  sub <- tempfile(fileext = ".csv"); map <- tempfile(fileext = ".csv")
  writeLines(c('"","Participant","genotype"', '"1","pt-a","T21"', '"2","pt-b","mosaic_T21"',
               '"3","pt-c","D21"'), sub)
  writeLines(c(",Participant ID,External Participant ID,HTP_num,Participant",
               "0,pt-a,HTP0001,0001,pt-a", "1,pt-b,HTP0002,0002,pt-b",
               "2,pt-c,KFDS999,,pt-c"), map)
  got <- karyotype_subtype_table(sub, map)
  expect_equal(nrow(got), 2)
  expect_equal(got[subject_id == "HTP0002", karyotype_subtype], "mosaic_T21")
  expect_setequal(names(got), c("subject_id", "karyotype_subtype"))
})

test_that("add_karyotype_subtype leaves unmapped samples NA and preserves row order and count", {
  meta <- data.table::data.table(LabID = c("HTP0001A", "HTP0002B2", "HTP0009A"),
                                 subject_id = c("HTP0001", "HTP0002", "HTP0009"))
  subs <- data.table::data.table(subject_id = c("HTP0002", "HTP0001"),
                                 karyotype_subtype = c("mosaic_T21", "T21"))
  got <- add_karyotype_subtype(meta, subs)
  expect_equal(nrow(got), 3)
  expect_equal(got$LabID, meta$LabID)
  expect_equal(got$karyotype_subtype, c("T21", "mosaic_T21", NA))
})
