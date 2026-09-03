# QC ABC-preEVADE ----------------------------------------------------
# library(OlinkAssay)
devtools::load_all()
# Background and Data Cleaning --------------------------------------------
# Background data from the received manifest
proj_dir <- here::here("tests/initial_data")

#TODO: remove after testing
npxData <- arrow::read_parquet(here::here(
  proj_dir,
  "FS19030710_Extended_NPX_2026-06-22.parquet"
))

npx_files <- c(
  here::here(proj_dir, "FS19030710_Extended_NPX_2026-06-22.parquet"),
  here::here(proj_dir, "FS19030718_Extended_NPX_2026-06-26.parquet"),
  here::here(proj_dir, "FS19030751_Extended_NPX_2026-06-22.parquet")
)

# read plates from disk
# because of the manifest, we need to specify `metadata_sheet` and `sample_column`
# As far as I know, the info in `extra_metadata` below has to be added manually here
# since it is not present in the files
# `project_column` would probably be find to be left to the default, but might as well
# sepcify it
objs <- purrr::map(
  .x = npx_files,
  .f = \(x) {
    OlinkAssayFromDisk(
      npx_file = x,
      metadata_file = here::here(
        proj_dir,
        "2026-05-15 James ABC-PreEVADE_Peds Serum_Olink (C.Guthridge)_Manifest.xlsx"
      ),
      metadata_sheet = "ManifestBuilder",
      sample_column = "Tube ID",
      project_column = "Project",
      extra_metadata = list(
        KitLot_info = "LC00358",
        PCLot_info = "1043509",
        SCLot_info = "1043510",
        NCLot_info = "1033075"
      ),
      verbose = FALSE
    )
  },
  .progress = TRUE
)

# combine the three plates
oa <- purrr::reduce(.x = objs, .f = concat)

# Performing batch correction using the median Plate Control Correction Procedure
oa <- batchCorrection(oa, method = "median")

# Performing Level 2 QC
oa <- scaleNPX(oa)
oa <- level2Process(oa)


# Example, categorical variable calls
rowData(oa) |>
  tibble::as_tibble(rownames = "Assay") |>
  dplyr::filter(assay_level_qc == "Categorical") |>
  dplyr::group_by(Assay) |>
  dplyr::summarise()

# Example, high variance variable calls
rowData(oa) |>
  tibble::as_tibble(rownames = "Assay") |>
  dplyr::filter(high_var_assay == "High Variance") |>
  dplyr::group_by(Assay) |>
  dplyr::summarise()

# Investigating Sample Failures
# abc_level_2 |>
#   dplyr::filter(SampleType == "SAMPLE") |>
#   dplyr::filter(is.na(ExtNPX)) |>
#   dplyr::group_by(Block, Assay) |>
#   dplyr::summarise() |>
#   dplyr::group_by(Block) |>
#   dplyr::summarise()

# Generation, Assay Performance Report
generate_performance_report(
  project_name = "ABC preEVADE Serum",
  correction_procedure = "Plate Control Batch Correction",
  sdp_directory = here::here("example", "SDP"),
  report_name = "ABC_preEVADE_assay_performance_report.pdf",
  level_1_file_path = stringr::str_glue("{proj_dir}/SDP/Level_1"),
  level_2_file_path = stringr::str_glue("{proj_dir}/SDP/Level_2")
)

# abc_level_2 |>
#   dplyr::group_by(Assay) |>
#   dplyr::summarise()

# Manifest Generation ----------------------------------------------------------------

# Preprocessing for manifest generation
abc_lvl1[["data"]] |>
  # Converting level 1 filename to the Assay_Filename column
  dplyr::bind_rows(.id = "Assay_Filename") |>
  # Filtering to only Technical Bridges and Samples
  dplyr::filter(SampleType == "SAMPLE") |>
  # Counting SampleQC entries per columnn
  dplyr::count(
    Assay_Filename,
    SampleID,
    Project,
    SampleQC
  ) |>
  # Entering SampleQC field Counts as separate columns
  tidyr::pivot_wider(
    names_from = SampleQC,
    values_from = n,
    values_fill = 0
  ) |>
  # Counting the percentage of Fail and Warn analytes per Sample
  dplyr::mutate(
    FAIL = dplyr::coalesce(FAIL, 0),
    WARN = dplyr::coalesce(WARN, 0),
    Analyte_Failure_Pct = 100 *
      FAIL /
      dplyr::n_distinct(abc_level_2[["Assay"]]),
    Analyte_Warn_Pct = 100 * WARN / dplyr::n_distinct(abc_level_2[["Assay"]])
  ) |>
  # Joining the manifest information from the original received manifest
  dplyr::left_join(
    manifest_information,
    by = dplyr::join_by(SampleID == `Tube ID`)
  ) |>
  # Fields not included in the original manifest
  dplyr::mutate(
    Same_Day_Processing = NA,
    Serum_Hemolysis = NA,
    Serum_Lipemia = NA,
    Processing_Comments = NA
  ) |>
  # Selecting and Ordering the Columns Prior to Writing the Manifest
  dplyr::select(
    SubjectRef,
    VisitRef,
    `Sample Alias`,
    `Sample Type`,
    Same_Day_Processing,
    Serum_Hemolysis,
    Serum_Lipemia,
    Processing_Comments,
    Assay_Filename,
    SampleID,
    Analyte_Failure_Pct,
    Analyte_Warn_Pct,
    Project
  ) |>
  # Reformatting Manifest Column Names
  dplyr::rename(
    Assay_Sample_ID = SampleID,
    Sample_Alias = `Sample Alias`,
    Sample_Type = `Sample Type`
  ) |>
  # Writing a preliminary excel template to be cleaned
  # according to the format of the manifest_template.xlsx file in the resources folder

  # TODO: need to make the output filename either an argument somewhere or
  # dependent on the projects included
  readr::write_csv("manifest_raw.csv")


# Investigating Potential Dispensation Effect ------------------------------------
# Example of potential challenge
plate_info <- readxl::read_excel(
  here::here(
    "FS19030710_NPX_2026-06-22.xlsx"
  ),
  sheet = "Sample Information"
)

failed <- abc_level_2 |>
  dplyr::filter(AssayType == "assay", AssayQC != "WARN") |>
  dplyr::filter(SampleQC == "FAIL") |>
  dplyr::pull(SampleID)

# Analytes with no entries among all samples
# Highlighting that the ExtNPX field is stricly from Manufacturer Calls
abc_level_2 |>
  dplyr::filter(is.na(ExtNPX)) |>
  dplyr::group_by(Assay) |>
  dplyr::summarise(count = n()) |>
  dplyr::filter(count != 11)

# Highlighting the samples with higher than expected analyte failure
failed <- abc_level_2 |>
  dplyr::filter(SampleType == "SAMPLE") |>
  dplyr::filter(Assay == "ABCF3") |>
  dplyr::filter(is.na(ExtNPX))

# Investigating the 96 and 384 Well Plate Positions of the sample failures
# To look at potential patterns in sample failure
plate_info |>
  dplyr::filter(sample_id %in% failed$SampleID)
