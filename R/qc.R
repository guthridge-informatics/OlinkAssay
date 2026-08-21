# TODO: I don't think any of these functions should have a side effect
# where they save a file. Since these all are saving the same data
# that they return, we just yank that and put it into a new file

#' olink_lvl1
#' Level 1 QC
#'
#' @param olink_files Output from [olink_reader()]. A list with two members:
#' 1. `"data"`: a named list of [`tibbles`][tibble::tibble] containing Olink data
#' 2. `"manifest"`: a named list of [`tibbles`][tibble::tibble] containing manifest data
#' @param proj_names Character vector with the projects the data is associated with
#' @param KitLot_info Lot number for the Olink Explore HT Kit used for the run
#' @param PCLot_info Lot number for Plate Control used for the Olink Explore HT run
#' @param SCLot_info Unknown
#' @param NCLot_info Lot number for Negative Control used for Olink Explore HT run
#' @param proj_dir path to where the files should be saved in the standard data package output format
#'
#' @returns
#'
#' @importFrom purrr map
#' @importFrom dplyr filter left_join join_by mutate
#' @importFrom stringr str_glue
#'
#' @export
#' @examples
olink_lvl1 <- function(
  olink_files,
  proj_names,
  KitLot_info,
  PCLot_info,
  SCLot_info,
  NCLot_info,
  proj_dir = NULL
) {
  proj_dir <- proj_dir %||% getwd()

  names <- names(olink_files[["data"]]) # get the names of all the files aka list names
  # store in the standard control names in the SampleID column
  ctrls <- c(
    "PC1",
    "PC2",
    "PC3",
    "PC4",
    "PC5",
    "NC1",
    "NC2",
    "SC1",
    "SC2",
    "SC3"
  )
  # as far as I can tell, there should only be one manifest at a time
  # so why it is a named list and we have to iterate through it seems pointless?
  # and can cause issues, e.g. trying to map a single entry results in a "can't do this
  # to a character vector" error
  manifest_filtered <-
    purrr::map(
      .x = olink_files[["manifest"]],
      .f = \(x) dplyr::filter(x, Project %in% proj_names)
    ) # filter the manifest file by the project names first
  multifile_write(
    .data = manifest_filtered,
    file_extension = "csv",
    proj_dir = proj_dir
  ) # write the SMDs

  # filter the current raw data by only selecting the sample only pertaining to the project names
  data_filtered <-
    purrr::map(
      .x = names,
      .f = \(x) {
        dplyr::filter(
          olink_files[["data"]][[x]],
          SampleID %in% c(ctrls, manifest_filtered[[1]][["SampleID"]])
        ) |>
          dplyr::left_join(
            olink_files[["manifest"]][[1]],
            by = dplyr::join_by("SampleID" == "SampleID")
          ) |>
          dplyr::mutate(
            KitLot = KitLot_info,
            PCLot = PCLot_info,
            SCLot = SCLot_info,
            NCLot = NCLot_info,
            AssayTechIssue = NA,
            UniqueID = stringr::str_glue("{PlateID}-{WellID}-{SampleID}")
          )
      }
    ) # filter only the sample IDs belongs to the projects
  names(data_filtered) <- names
  multifile_write(
    .data = data_filtered,
    file_extension = "parquet",
    proj_dir = proj_dir
  ) # write the Level 1 parquet

  list(
    data = data_filtered,
    manifest = manifest_filtered
  )
}


#' @title olink_qc
#' @description tallies the number of failed and passed assay and samples per each 96-well plate
#'
#' @param data [`tibble`][`tibble::tibble`] containing Olink data
#'
#' @returns Named list with
#' 1. `AssayOlinkQC`: table with assay-level quality control counts
#' 2. `SampleOlinkQC`: table with sample-level quality control counts
#'
#' @importFrom dplyr filter group_by tally ungroup
#' @export
#' @examples
olink_qc <- function(.data) {
  table_assay <-
    .data |>
    dplyr::filter(AssayType == "assay") |> # filter for only assays instead of all extension, plate controls
    dplyr::group_by(PlateID, OlinkID, AssayQC) |> # group by PlateID and OlinkID to get the split among the Assay QC
    dplyr::tally() |> # tabulate the counts
    dplyr::ungroup() |> # this is to allow the ungrouping of the tibble and return to orignal tibble without grouping
    dplyr::group_by(PlateID, AssayQC) |> # grouping to allow tabulation of assays
    dplyr::tally()

  table_sample <-
    .data |>
    dplyr::filter(SampleType == "SAMPLE") |> # filter for only assays instead of all extension, plate controls
    dplyr::group_by(PlateID, SampleID, SampleQC) |> # group by PlateID and OlinkID to get the split among the Assay QC
    dplyr::tally() |> # tabulate the counts
    dplyr::mutate(Frequency = prop.table(n)) # calculate the frequencies within each sample

  list(
    AssayOlinkQC = table_assay,
    SampleOlinkQC = table_sample
  )
}

## Level 2 QC Part 1
#' @title olink_lvl2_prep
#' @description
#'
#' @param data a [`tibble`][`tibble::tibble`] containing batch-corrected Olink data
#'
#' @returns [`tibble::tibble`]
#'
#' @importFrom dplyr mutate filter group_by summarise case_when select left_join
#' @export
#' @examples
olink_lvl2_prep <- function(.data) {
  .data <-
    .data |>
    dplyr::mutate(
      LogProtExp = ExtNPX_Corrected + log2(1e5), # transform into LogProExp
      LogProtExp_Raw = ExtNPX_Corrected + log2(1e5) # transform into LogProExp
    )

  ht_nc_vals <- .data |>
    dplyr::filter(SampleType == "NEGATIVE_CONTROL") |>
    dplyr::group_by(Assay, OlinkID) |>
    dplyr::summarise(
      median_nc = median(na.omit(LogProtExp)),
      iqr_nc = as.numeric(quantile(na.omit(LogProtExp), 0.75)),
      .groups = 'drop'
    )

  # Calculating Plate Control coefficient of variance
  ht_pc_vals <- .data |>
    dplyr::filter(SampleType == "PLATE_CONTROL") |>
    dplyr::group_by(Assay, OlinkID) |>
    dplyr::summarise(
      pc_cv = 100 * sd(LogProtExp) / mean(LogProtExp),
      .groups = 'drop'
    ) |>
    dplyr::mutate(
      high_var_assay = dplyr::case_when(
        pc_cv > 20 ~ "High Variance",
        TRUE ~ "Pass"
      )
    ) |>
    dplyr::select(-pc_cv)

  # This is the "sample level" qc, calculates ith sample in jth assay that needs to be replaced with zero or LLOQ
  # also labels those values in a new column - sample_level_qc
  .data |>
    dplyr::filter(SampleType == "SAMPLE") |>
    dplyr::left_join(ht_nc_vals, by = c("Assay", "OlinkID")) |>
    dplyr::left_join(ht_pc_vals, by = c("Assay", "OlinkID")) |>
    dplyr::mutate(
      sample_level_qc = dplyr::case_when(
        LogProtExp < median_nc ~ "Below LLOD",
        LogProtExp < iqr_nc ~ "Below LLOQ",
        TRUE ~ "Pass"
      )
    ) |>
    dplyr::mutate(
      LogProtExp = dplyr::case_when(
        LogProtExp < median_nc ~ 0,
        LogProtExp < iqr_nc ~ iqr_nc,
        TRUE ~ LogProtExp
      )
    )
}

#' @title olink_lvl2
#' @description Level 2 QC Part 2
#'
#' @param data tibble or list of tibbles
#' @param multifile bool
#' @param proj_dir path to save to
#'
#' @returns
#'
#' @importFrom dplyr group_by summarise mutate case_when select left_join n_distinct n
#' @importFrom tidyr pivot_wider
#' @importFrom arrow write_parquet
#' @importFrom purrr list_rbind
#'
#' @export
#' @examples
olink_lvl2 <- function(.data, multifile = TRUE, proj_dir = NULL) {
  proj_dir <- proj_dir %||% getwd()

  # concatenate the data
  if (is.list(.data)) {
    .data <- purrr::list_rbind(.data)
  }

  # specifying number of samples present in the total combined dataset
  n_samples <- .data |> dplyr::select(SampleID) |> dplyr::n_distinct()

  # Assay level QC - if 50% of samples are below LLOQ, labeled as semi-continuous
  # if 75% of samples are below LLOD, labeled as categorical
  # Test to adjust how I calculate categorical, semi-continuous, or continuous
  ht_scaled_npx_assay <- .data |>
    dplyr::group_by(Assay, OlinkID, sample_level_qc) |>
    dplyr::summarise(
      percentage = 100 * dplyr::n() / n_samples,
      .groups = 'drop'
    ) |>
    tidyr::pivot_wider(
      names_from = sample_level_qc,
      values_from = percentage
    ) |>
    dplyr::mutate(
      `Below LLOD` = dplyr::case_when(
        is.na(`Below LLOD`) ~ 0,
        TRUE ~ `Below LLOD`
      ),
      `Below LLOQ` = dplyr::case_when(
        is.na(`Below LLOQ`) ~ 0,
        TRUE ~ `Below LLOQ`
      )
    ) |>
    dplyr::mutate(
      `Below LLOQ` = `Below LLOQ` + `Below LLOD`,
      assay_level_qc = dplyr::case_when(
        `Below LLOD` > 75 ~ "Categorical",
        `Below LLOQ` > 50 ~ "Semi-Continuous",
        TRUE ~ "Continuous"
      )
    ) |>
    dplyr::select(-c(`Below LLOD`, `Below LLOQ`))

  data <- dplyr::left_join(
    .data,
    dplyr::select(ht_scaled_npx_assay, OlinkID, assay_level_qc),
    by = dplyr::join_by(OlinkID)
  )
  lvl2_output_dir <- stringr::str_glue("{proj_dir}/SDP/Level_2/")
  if (!lvl2_output_dir %in% list.dirs(proj_dir, recursive = TRUE)) {
    dir.create(path = lvl2_output_dir, recursive = TRUE)
  }
  arrow::write_parquet(
    data,
    sink = stringr::str_glue("{lvl2_output_dir}/Level_2_SDP.parquet")
  ) # write the parquet files of the level 2
  data
}
