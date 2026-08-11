testthat::test_that(desc = "ingest_olink_data", code = {
  expected_data <- readRDS(here::here(
    "tests",
    "testthat",
    "data",
    "common",
    "ingest_olink_data_expected.RDS"
  ))
  test_data_directory <- here::here("tests", "initial_data")
  testthat::expect_mapequal(
    olinkqc::ingest_olink_data(directory = test_data_directory),
    expected_data
  )
})

testthat::test_that(desc = "ingest_manifest", code = {
  expected_manifest <- readRDS(here::here(
    "tests",
    "testthat",
    "data",
    "io",
    "ingest_manifest_expected.RDS"
  ))
  test_data_directory <- here::here("tests", "initial_data")
  testthat::expect_mapequal(
    object = olinkqc::ingest_manifest(
      directory = test_data_directory,
      sample_column = "Tube ID",
      project_column = "Project",
      manifest_sheet = "ManifestBuilder"
    ),
    expected = expected_manifest
  )
})

testthat::test_that(desc = "multifile_write", code = {
  expected_data <- readRDS(here::here(
    "tests",
    "testthat",
    "data",
    "common",
    "ingest_olink_data_expected.RDS"
  ))
  expected_manifest <- readRDS(here::here(
    "tests",
    "testthat",
    "data",
    "io",
    "ingest_manifest_expected.RDS"
  ))
  path <- withr::local_tempfile()
  multifile_write(expected_data, file_extension = "parquet", proj_dir = path)
  multifile_write(expected_manifest, file_extension = "csv", proj_dir = path)
  testthat::expect_snapshot_file(
    stringr::str_glue(
      "{path}/SDP/Level_1/FS19030710_NPX_2026-06-22.parquet"
    )
  )
  testthat::expect_snapshot_file(
    stringr::str_glue(
      "{path}/SDP/Level_1/FS19030718_NPX_2026-06-26.parquet"
    )
  )
  testthat::expect_snapshot_file(
    stringr::str_glue(
      "{path}/SDP/Level_1/FS19030751_NPX_2026-06-22.parquet"
    )
  )
  testthat::expect_snapshot_file(
    stringr::str_glue(
      "{path}/SDP/Level_1/2026-05-15 James ABC-PreEVADE_Peds Serum_Olink (C.Guthridge)_Manifest.csv"
    )
  )
})
