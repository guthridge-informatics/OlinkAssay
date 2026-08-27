#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @param .data stuff
#' @param ... other arguments to pass on
.OlinkAssay <- setClass(
  Class = "OlinkAssay",
  contains = "SummarizedExperiment"
)

#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
OlinkAssay <- function(
  npxData,
  colData = NULL,
  assay_tables = c(
    "Count",
    "ExtNPX",
    "PCNormalizedNPX",
    "NPX"
  ),
  ...
) {
  assay_medians <- ctrl_ref(npxData)
  global_medians <- global_ref(npxData)

  .metadata <- dplyr::select(
    npxData,
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

  .metadata[["AssayMedians"]] <- assay_medians
  .metadata[["GlobalMedians"]] <- global_medians

  .metadata[["SampleBlockQCWarn"]] <-
    dplyr::filter(
      .data = npxData,
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
      .data = npxData,
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
        npxData |>
          dplyr::select(
            -Panel,
            -SoftwareVersion,
            -SoftwareName,
            -PanelDataArchiveVersion,
            -PreProcessingVersion,
            -PreProcessingSoftware,
            -InstrumentType
          ) |>
          .extract_dfs(
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
    npxData |>
    .extract_dfs(
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
    dplyr::distinct()

  if (!is.null(colData)) {
    dplyr::left_join(
      x = .coldata,
      y = colData,
      by = dplyr::join_by("SampleID" == "SampleID")
    )
  }
  .coldata <-
    tibble::column_to_rownames(.data = .coldata, var = "SampleID") |>
    S4Vectors::DataFrame()

  .rowdata <-
    npxData |>
    .extract_dfs(
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

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = sample_df,
    rowData = .rowdata,
    colData = .coldata,
    metadata = .metadata
  )
  .OlinkAssay(se)
}

S4Vectors::setValidity2("OlinkAssay", function(object) {
  msg <- NULL

  required_assays <- c("Count", "ExtNPX", "NPX")
  for (x in required_assays) {
    if (!x %in% SummarizedExperiment::assayNames(object)) {
      msg <- c(msg, glue::glue("{x} is a required assay"))
    }
  }

  if (is.null(msg)) {
    TRUE
  } else {
    msg
  }
})

#' @title assay
#' @description Get assay data for OlinkAssay object
#'
#' @param i name of the assay in `assays(x)` to retreive. If no name is provided
#'    attempts to extract the assay named (in decending order) "ExtNPX_Corrected",
#'.   "ExtNPX", "PCNormalizedNPX", "Correction", "NPX", "Count"
#' @param withDimnames retain row and column names? (Default: TRUE)
#'
#' @inheritParams SummarizedExperiment::assay
#' @export
setGeneric("assay", function(x, ...) standardGeneric("assay"))

#' @export
#' @importFrom SummarizedExperiment assay
setMethod(
  f = "assay",
  signature = "OlinkAssay",
  definition = function(
    x,
    i = NULL,
    withDimnames = TRUE
  ) {
    if (is.null(i)) {
      i <- intersect(
        c(
          "ExtNPX",
          "PCNormalizedNPX",
          "NPX",
          "Count"
        ),
        colnames(rowData(x))
      )[[1]]
    }
    SummarizedExperiment::assay(x = x, i = i, withDimnames = withDimnames)
  }
)

#' @export
setGeneric("npx", function(x, ...) standardGeneric("npx"))

#' @export
#' @importFrom SummarizedExperiment assay
setMethod(
  f = "npx",
  signature = "OlinkAssay",
  definition = function(x, withDimnames = TRUE) {
    SummarizedExperiment::assay(x = x, i = "NPX", withDimnames = withDimnames)
  }
)


#' @export
setGeneric("colData", function(x, ...) standardGeneric("colData"))

#' @export
#' @importFrom methods slot
setMethod(
  f = "colData",
  signature = "OlinkAssay",
  definition = function(
    x,
    withDimnames = TRUE
  ) {
    slot(x, "colData")
  }
)

#' @export
setGeneric("rowData", function(x, ...) standardGeneric("rowData"))

#' @export
#' @importFrom S4Vectors mcols
setMethod(
  f = "rowData",
  signature = "OlinkAssay",
  definition = function(
    x,
    withDimnames = TRUE
  ) {
    S4Vectors::mcols(x)
  }
)


#' @export
#' @importFrom arrow read_parquet
readFromDisk <- function(
  npx_file = NULL,
  metadata_file = NULL,
  metadata_sheet = 1,
  sample_column = "SampleID",
  project_column = "Project",
  additional_columns = NULL
) {
  if (is.character(npx_file)) {
    df <- arrow::read_parquet(npx_file)
  } else {
    df <- S4Vectors::DataFrame()
  }

  if (is.character(metadata_file)) {
    ext <- stringr::str_split_i(
      string = metadata_file,
      pattern = .Platform$file.sep,
      i = -1
    ) |>
      stringr::str_split_i(pattern = "\\.", i = -1)

    if (!is.null(sample_column)) {
      sample_column <- rlang::sym(sample_column)
    }
    if (!is.null(project_column)) {
      project_column <- rlang::sym(project_column)
    }

    if (!is.null(additional_columns)) {
      loc <- tidyselect::eval_select(
        rlang::expr(additional_columns),
        data = manifest
      )
    } else {
      loc <- NULL
    }

    md <- switch(
      EXPR = ext,
      "xlsx" = readxl::read_excel(metadata_file, sheet = metadata_sheet),
      "csv" = readr::read_csv(metadata_file)
    ) |>
      dplyr::select(
        SampleID = {{ sample_column }},
        Project = {{ project_column }},
        loc
      )
  } else {
    md <- NULL
  }

  se <- OlinkAssay(npxData = df, colData = md)
}

#' @export
setGeneric("medianCorrection", function(x, ...) {
  standardGeneric("medianCorrection")
})

#' @export
setMethod(
  f = "medianCorrection",
  signature = "OlinkAssay",
  definition = function(
    x,
    method = c("median", "global median")
  ) {
    method <- match.arg(method)

    correction_values <- switch(
      EXPR = method,
      "median" = metadata(x)[["AssayMedians"]],
      "global median" = metadata(x)[["GlobalMedians"]]
    )

    median_means <-
      correction_values |>
      dplyr::select(-Variance) |>
      tidyr::pivot_wider(
        names_from = "PlateRef",
        values_from = "Median"
      ) |>
      dplyr::rowwise() |>
      dplyr::transmute(
        OlinkID = OlinkID,
        ReferenceMedian = mean(dplyr::c_across(
          -tidyselect::contains("OlinkID")
        )),
      )

    meds_correction <-
      dplyr::left_join(
        x = correction_values,
        y = median_means,
        by = dplyr::join_by(OlinkID)
      ) |>
      dplyr::mutate(Correction = Median - ReferenceMedian) |>
      dplyr::select(-Median, -Variance)

    olinkid_to_assay <- rowData(x) |>
      tibble::as_tibble(rownames = "assay") |>
      dplyr::select(OlinkID, assay)

    sample_to_plate <- colData(x) |>
      tibble::as_tibble(rownames = "SampleID") |>
      dplyr::mutate(
        PlateRef = stringr::str_remove(string = PlateID, pattern = "_.*$")
      ) |>
      dplyr::select(SampleID, PlateRef)

    npx_values <- assay(x, i = "ExtNPX") |>
      tibble::rownames_to_column(var = "assay") |>
      tidyr::pivot_longer(
        -assay,
        names_to = "SampleID",
        values_to = "ExtNPX"
      )

    data_correction <-
      dplyr::left_join(
        npx_values,
        sample_to_plate,
        by = dplyr::join_by(SampleID)
      ) |>
      dplyr::left_join(
        y = olinkid_to_assay,
        by = dplyr::join_by(assay)
      ) |>
      dplyr::left_join(
        y = meds_correction,
        by = dplyr::join_by(OlinkID, PlateRef)
      ) |>
      dplyr::mutate(ExtNPX_Corrected = ExtNPX - Correction)

    assay(x, i = "ExtNPX_corrected") <- data_correction |>
      dplyr::select(SampleID, assay, ExtNPX_Corrected) |>
      dplyr::distinct() |>
      tidyr::pivot_wider(
        names_from = "SampleID",
        values_from = "ExtNPX_Corrected"
      ) |>
      tibble::column_to_rownames("assay")

    metadata(x)[["BatchCorrectionMethod"]] <- method

    x
  }
)
