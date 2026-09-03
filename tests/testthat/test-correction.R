testthat::test_that(desc = "ctrl_ref", code = {
  example_data <- readRDS(here::here(
    "tests",
    "testthat",
    "data",
    "correction",
    "example_single_olink_data_tbl.RDS"
  ))
  expected_correction <- readRDS(here::here(
    "tests",
    "testthat",
    "data",
    "correction",
    "expected_ctrl_ref_tbl.RDS"
  ))
  testthat::expect_equal(
    olinkqc::ctrl_ref(.data = example_data),
    expected_correction
  )
})

testthat::test_that(desc = "global_ref", code = {
  example_data <- readRDS(here::here(
    "tests",
    "testthat",
    "data",
    "correction",
    "example_single_olink_data_tbl.RDS"
  ))
  expected_correction <- readRDS(here::here(
    "tests",
    "testthat",
    "data",
    "correction",
    "expected_global_ref_tbl.RDS"
  ))
  testthat::expect_equal(
    olinkqc::global_ref(.data = example_data),
    expected_correction
  )
})

# I don't think we need a particular test for `median_correction` since the one below
# pretty much covers it
testthat::test_that(desc = "batch_correction", code = {
  example_data <- readRDS(here::here(
    "tests",
    "testthat",
    "data",
    "common",
    "expected_ingest_olink_data.RDS"
  ))
  expected_median_corrected <- readRDS(here::here(
    "tests",
    "testthat",
    "data",
    "correction",
    "expected_median_corrected.RDS"
  ))
  expected_global_median_corrected <- readRDS(here::here(
    "tests",
    "testthat",
    "data",
    "correction",
    "expected_global_median_corrected.RDS"
  ))

  testthat::expect_equal(
    object = batch_correction(.data = example_data, method = "median"),
    expected = expected_median_corrected
  )
  testthat::expect_equal(
    object = batch_correction(.data = example_data, method = "global"),
    expected = expected_global_median_corrected
  )
})
