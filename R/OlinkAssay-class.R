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

OlinkAssay <- function(
  assays = SimpleList(),
  rowData = NULL,
  rowRanges = NULL,
  colData = DataFrame(),
  metadata = list(),
  checkDimnames = TRUE
) {
  .OlinkAssay(
    SummarizedExperiment(
      assays = assays,
      rowData = rowData,
      rowRanges = rowRanges,
      colData = colData,
      metadata = metadata,
      checkDimnames = checkDimnames
    )
  )
}

#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
OlinkAssayFromNPX <- function(
  npxData,
  colData = NULL,
  assay_tables = c(
    "Count",
    "ExtNPX",
    "PCNormalizedNPX",
    "NPX"
  ),
  verbose = FALSE,
  ...
) {
  if (verbose) {
    message("Calculating assay medians...")
  }
  assay_medians <- S4Vectors::DataFrame(ctrl_ref(npxData))
  if (verbose) {
    message("Calculating global medians...")
  }
  global_medians <- S4Vectors::DataFrame(global_ref(npxData))

  if (verbose) {
    message("Extracting metadata")
  }
  plate_ref <- stringr::str_remove(
    string = npxData[["PlateID"]],
    pattern = "_.*$"
  ) |>
    unique()

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
    as.list() |>
    # look, if you've got a better way to name a list using a variable
    # but actually evaluate the variable I'd love to see it
    purrr::map(.f = \(x) {
      item <- list(x)
      names(item) <- plate_ref
      item
    })
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

  if (verbose) {
    message("Extracting assay data")
  }
  sample_df <-
    purrr::map(
      .x = assay_tables,
      .f = \(x) {
        if (verbose) {
          message(glue::glue("Pivoting {x}"))
        }
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

  if (verbose) {
    message("Extracting coldata")
  }
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
    if (verbose) {
      message("Adding external metadata")
    }
    dplyr::left_join(
      x = .coldata,
      y = colData,
      by = dplyr::join_by("SampleID" == "SampleID")
    )
  }
  .coldata <-
    tibble::column_to_rownames(.data = .coldata, var = "SampleID") |>
    S4Vectors::DataFrame()

  if (verbose) {
    message("Extracting rowdata")
  }
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

  if (verbose) {
    message("Creating object")
  }
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
      tibble::as_tibble() |>
      dplyr::select(-Variance) |>
      tidyr::pivot_wider(
        names_from = "PlateRef",
        names_glue = "{PlateRef}_median",
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
        x = tibble::as_tibble(correction_values),
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

    npx_values <- SummarizedExperiment::assay(x, i = "ExtNPX") |>
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

#' @export
setGeneric("concat", function(x, ...) standardGeneric("concat"))

#' @export
#' @importFrom S4Vectors mcols
#' @importFrom stats setNames
setMethod(
  f = "concat",
  signature = "OlinkAssay",
  definition = function(
    x,
    y,
    batch_correction_method = NULL
  ) {
    concat_coldata <- rbind(colData(x), colData(y))
    if (all(rownames(rowData(x)) == rownames(rowData(y)))) {
      concat_rowdata <- rowData(x)
    } else {
      missing_in_x <- rownames(rowData(y))[
        !rownames(rowData(y)) %in% rownames(rowData(x))
      ]
      concat_rowdata <- rbind(rowData(x), rowData(y)[, missing_in_x])
    }

    # merge metadata
    single_value_cols <- c(
      "Panel",
      "SoftwareVersion",
      "SoftwareVersion",
      "SoftwareName",
      "PanelDataArchiveVersion",
      "PreProcessingVersion",
      "PreProcessingSoftware",
      "InstrumentType"
    )
    md_list <- purrr::map(
      .x = single_value_cols,
      .f = \(z) {
        c(metadata(x)[[z]], metadata(y)[[z]])
      }
    ) |>
      stats::setNames(single_value_cols)

    md_list[["AssayMedians"]] <- rbind(
      metadata(x)[["AssayMedians"]],
      metadata(y)[["AssayMedians"]]
    )
    md_list[["GlobalMedians"]] <- rbind(
      metadata(x)[["GlobalMedians"]],
      metadata(y)[["GlobalMedians"]]
    )

    # `base_assays` being those that are present when naive data
    # is first imported. If we are merging data, I am not yet sure
    # that merging already corrected data should be included and am
    # instead relying on repeating any previous corrections.
    # will have to test
    base_assays <- c("Count", "ExtNPX", "PCNormalizedNPX", "NPX")
    merged_assays <- purrr::map(
      .x = base_assays,
      .f = \(z) {
        merge(assay(x, i = z), assay(y, i = z), by = "row.names") |>
          tibble::column_to_rownames("Row.names")
      }
    ) |>
      purrr::set_names(base_assays)

    concat_oa <- OlinkAssay(
      assays = merged_assays,
      colData = concat_coldata,
      rowData = `rownames<-`(concat_rowdata, NULL),
      metadata = md_list
    )

    if (!is.null(batch_correction_method)) {
      concat_oa <- medianCorrection(concat_oa, method = batch_correction_method)
    }

    concat_oa
  }
)

reshape_medians_tbl <- function(.data) {
  .data |>
    tibble::as_tibble() |>
    dplyr::select(-Variance) |>
    tidyr::pivot_wider(
      names_from = "PlateRef",
      names_glue = "{PlateRef}_median",
      values_from = "Median"
    )
}
