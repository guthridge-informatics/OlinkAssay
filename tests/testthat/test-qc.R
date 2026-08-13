testthat::test_that(desc = "Olink_lvl1", code = {
  example_olink_reader_results <- readRDS(
    here::here(
      "tests",
      "testthat",
      "data",
      "common",
      "example_olink_reader_results.RDS"
    )
  )
  expected_olink_lvl1_output <- readRDS(
    here::here(
      "tests",
      "testthat",
      "data",
      "qc",
      "expected_olink_lvl1_results.RDS"
    )
  )
  tmp_path <- withr::local_tempfile()
  lvl1_results <-
    olinkqc::Olink_lvl1(
      olink_files = example_olink_reader_results,
      proj_names = c("Bridge", "ABC-preEVADE"),
      KitLot_info = "LC00358",
      PCLot_info = "1043509",
      SCLot_info = "1043510",
      NCLot_info = "1033075",
      proj_dir = tmp_path
    )
  testthat::expect_equal(
    object = lvl1_results,
    expected = expected_olink_lvl1_output
  )
  testthat::expect_snapshot_file(
    stringr::str_glue(
      "{tmp_path}/SDP/Level_1/FS19030710_NPX_2026-06-22.parquet"
    )
  )
  testthat::expect_snapshot_file(
    stringr::str_glue(
      "{tmp_path}/SDP/Level_1/FS19030718_NPX_2026-06-26.parquet"
    )
  )
  testthat::expect_snapshot_file(
    stringr::str_glue(
      "{tmp_path}/SDP/Level_1/FS19030751_NPX_2026-06-22.parquet"
    )
  )
  testthat::expect_snapshot_file(
    stringr::str_glue(
      "{tmp_path}/SDP/Level_1/2026-05-15 James ABC-PreEVADE_Peds Serum_Olink (C.Guthridge)_Manifest.csv"
    )
  )
})

testthat::test_that(desc = "olink_qc", code = {
  example_olink_data_tbl <- readRDS(here::here(
    "tests",
    "testthat",
    "data",
    "correction",
    "example_single_olink_data_tbl.RDS"
  ))
  expected_olink_qc_output <- readRDS(
    here::here(
      "tests",
      "testthat",
      "data",
      "qc",
      "expected_olink_qc_results.RDS"
    )
  )
  olink_qc_results <- Olink_qc(.data = example_olink_data_tbl)
  testthat::expect_equal(
    object = olink_qc_results,
    expected = expected_olink_qc_output
  )
})

testthat::test_that(desc = "Olink_lvl2_prep", code = {
  example_batch_corrected_olink_data_tbl <- readRDS(
    here::here(
      "tests",
      "testthat",
      "data",
      "qc",
      "example_batch_corrected_olink_data_tbl.RDS"
    )
  )

  expected_olink_lvl2_prep <- readRDS(
    here::here(
      "tests",
      "testthat",
      "data",
      "qc",
      "expected_olink_lvl2_prep_results.RDS"
    )
  )
  testthat::expect_equal(
    object = Olink_lvl2_prep(example_batch_corrected_olink_data_tbl),
    expected = expected_olink_lvl2_prep
  )
})

testthat::test_that(desc = "Olink_lvl2", code = {
  example_lvl2_prepped_results <- readRDS(
    here::here(
      "tests",
      "testthat",
      "data",
      "qc",
      "example_olink_lvl2_prep_results.RDS"
    )
  )
  expected_olink_lvl2_results <- readRDS(
    here::here(
      "tests",
      "testthat",
      "data",
      "qc",
      "expected_olink_lvl2_results.RDS"
    )
  )
  tmp_path <- withr::local_tempfile()
  lvl2_results <- olinkqc::Olink_lvl2(
    .data = example_lvl2_prepped_results,
    multifile = TRUE,
    proj_dir = tmp_path
  )
  testthat::expect_equal(
    object = lvl2_results,
    expected = expected_olink_lvl2_results
  )
  testthat::expect_snapshot_file(
    stringr::str_glue(
      "{tmp_path}/SDP/Level_2/Level_2_SDP.parquet"
    )
  )
})
