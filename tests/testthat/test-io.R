testthat::test_that(desc = "ingest_olink_data", code = {
  expected_ingest_olink_data <- readRDS(here::here(
    "tests",
    "testthat",
    "data",
    "common",
    "expected_ingest_olink_data.RDS"
  ))
  test_data_directory <- here::here("tests", "initial_data")
  testthat::expect_mapequal(
    object = olinkqc::import_olink_data(input = test_data_directory),
    expected = expected_ingest_olink_data
  )
})

testthat::test_that(desc = "ingest_manifest", code = {
  expected_ingest_manifest <- readRDS(here::here(
    "tests",
    "testthat",
    "data",
    "io",
    "expected_ingest_manifest.RDS"
  ))
  test_data_directory <- here::here("tests", "initial_data")
  testthat::expect_mapequal(
    object = olinkqc::import_manifest(
      input = test_data_directory,
      sample_column = "Tube ID",
      project_column = "Project",
      manifest_sheet = "ManifestBuilder"
    ),
    expected = expected_ingest_manifest
  )
})

# TODO: use withr (I guess?) to clean up the written files
# even though they are in a temp dir, it seems like they persist
testthat::test_that(desc = "multifile_write", code = {
  expected_data <- readRDS(here::here(
    "tests",
    "testthat",
    "data",
    "common",
    "expected_ingest_olink_data.RDS"
  ))
  expected_manifest <- readRDS(here::here(
    "tests",
    "testthat",
    "data",
    "io",
    "expected_ingest_manifest.RDS"
  ))
  tmp_path <- withr::local_tempfile()
  multifile_write(
    .data = expected_data,
    file_extension = "parquet",
    proj_dir = tmp_path
  )
  multifile_write(
    .data = expected_manifest,
    file_extension = "csv",
    proj_dir = tmp_path
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
