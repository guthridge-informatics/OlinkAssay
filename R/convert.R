#' @title as_SummarizedExperiment
#' @description Convert the table from an Olink-HT parquet file into a
#'  SummarizedExperiment object. The data must already read from disk and be in memory.
#' @param .data [`tibble::tibble`] of raw Olink-HT data, directly from `arrow::read_parquet`
#' @param assay_tables Which columns contain data that should be retained as a feature-by-sample matrix in the `assays`
#'  slot (default: `Count`, `ExtNPX`, `ExtNPX_Corrected`, `Correction`, `NPX`, `PCNormalizedNPX`)
#' @param correction_method Which method should be used to correct batches, `median`, `global`, or `none` (default: `median`)?
#'
#' @returns A SummarisedExperiment version of the Olink data
#'
#' @export
#' @examples

as_SummarizedExperiment <- function(
  .data,
  assay_tables = c(
    "Count",
    "ExtNPX",
    "ExtNPX_Corrected",
    "Correction",
    "NPX",
    "PCNormalizedNPX"
  ),
  correction_method = c("median", "global", "none")
) {
  correction_method <- match.arg(correction_method)
  if (correction_method %in% c("median", "global")) {
    .data <- batch_correction.tbl_df(.data, method = correction_method)
  } else {
    assay_tables <- assay_tables[stringr::str_detect(
      string = assay_tables,
      pattern = "Correct",
      negate = TRUE
    )]
  }

  .metadata <- dplyr::select(
    .data,
    Panel,
    SoftwareVersion,
    SoftwareName,
    PanelDataArchiveVersion,
    PreProcessingVersion,
    PreProcessingSoftware,
    InstrumentType
  ) |>
    dplyr::distinct() |>
    tidyr::pivot_longer(
      cols = tidyselect::everything()
    ) |>
    tibble::deframe() |>
    as.list()

  .metadata[["SampleBlockQCWarn"]] <-
    dplyr::filter(
      .data = .data,
      SampleType == "SAMPLE",
      AssayType == "assay"
    ) |>
    dplyr::select(
      SampleID,
      Assay,
      SampleBlockQCWarn
    ) |>
    tidyr::pivot_wider(
      names_from = "SampleID",
      values_from = "SampleBlockQCWarn"
    ) |>
    tibble::column_to_rownames(var = "Assay")

  .metadata[["SampleBlockQCFail"]] <-
    dplyr::filter(
      .data = .data,
      SampleType == "SAMPLE",
      AssayType == "assay"
    ) |>
    dplyr::select(
      SampleID,
      Assay,
      SampleBlockQCFail
    ) |>
    tidyr::pivot_wider(
      names_from = "SampleID",
      values_from = "SampleBlockQCFail"
    ) |>
    tibble::column_to_rownames(var = "Assay")

  sample_df <-
    purrr::map(
      .x = assay_tables,
      .f = \(x) {
        .data |>
          dplyr::select(
            -Panel,
            -SoftwareVersion,
            -SoftwareName,
            -PanelDataArchiveVersion,
            -PreProcessingVersion,
            -PreProcessingSoftware,
            -InstrumentType
          ) |>
          extract_dfs(
            filter_assay_controls = TRUE,
            filter_sample_controls = TRUE
          ) |>
          purrr::pluck("SAMPLE") |>
          dplyr::select(SampleID, Assay, tidyselect::all_of(x)) |>
          tidyr::pivot_wider(names_from = "SampleID", values_from = x) |>
          tibble::column_to_rownames(var = "Assay")
      }
    ) |>
    purrr::set_names(assay_tables)

  .coldata <-
    .data |>
    extract_dfs(
      filter_assay_controls = TRUE,
      filter_sample_controls = TRUE
    ) |>
    purrr::pluck("SAMPLE") |>
    dplyr::select(
      SampleID,
      WellID,
      PlateID,
      OSICategory,
      OSISummary,
      OSITimeToCentrifugation,
      OSIPreparationTemperature
    ) |>
    dplyr::distinct() |>
    tibble::column_to_rownames(var = "SampleID") |>
    S4Vectors::DataFrame()

  .rowdata <-
    .data |>
    extract_dfs(
      filter_assay_controls = TRUE,
      filter_sample_controls = TRUE
    ) |>
    purrr::pluck("SAMPLE") |>
    dplyr::select(
      UniProt,
      Assay,
      AssayType,
      OlinkID,
      Normalization
      # AssayQC,
      # AssayQCWarn
    ) |>
    dplyr::distinct() |>
    tibble::column_to_rownames(var = "Assay") |>
    S4Vectors::DataFrame()

  SummarizedExperiment::SummarizedExperiment(
    assays = sample_df,
    rowData = .rowdata,
    colData = .coldata,
    metadata = .metadata
  )
}
