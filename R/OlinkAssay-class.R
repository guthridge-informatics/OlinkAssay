#### class definition ####
#' @title OlinkAssay
#' @description An S4 class to represent data from one or more Olink HT assays
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @slot assays a [`SimpleList`][`S4Vectors::SimpleList`] of feature-by sample
#'  [`DataFrames`][`S4Vectors::DataFrame`]. Must have
#'  "Count", "ExtNPX", "NPX" data
#' @slot rowData [`DataFrame`][`S4Vectors::DataFrame`] with feature metadata
#' @slot colData [`DataFrame`][`S4Vectors::DataFrame`] with sample metadata
#' @slot metadata list with assay-level metadata. Data is retained per-plate/run.
#'  Per-plate assay means are stored here
#'  under "AssayMedians" and "GlobalMedians"
#'
#' @name OlinkAssay
#' @export
#'
setClass(
  Class = "OlinkAssay",
  contains = "SummarizedExperiment",
  representation(
    negativeControls = "DataFrame",
    plateControls = "DataFrame"
  ),
  prototype(
    SummarizedExperiment(),
    negativeControls = new("DFrame"),
    plateControls = new("DFrame")
  )
)

#### constructors ####
new_OlinkAssay <- function(
  names,
  assays,
  rowData,
  colData,
  metadata,
  negativeControls,
  plateControls
) {
  new(
    "OlinkAssay",
    NAMES = names,
    assays = assays,
    elementMetadata = rowData,
    colData = colData,
    metadata = metadata,
    negativeControls = negativeControls,
    plateControls = plateControls
  )
}

OlinkAssay <- function(
  assays = Assays(),
  rowData = DataFrame(),
  colData = DataFrame(),
  metadata = list(),
  negativeControls = DataFrame(),
  plateControls = DataFrame()
) {
  first_assay <- assays@data[[1]]
  if (!all(colnames(first_assay) == rownames(colData))) {
    stop(
      "Names (likely SampleIDs) do not match between the Assay columns and colData rows"
    )
  }
  if (!all(rownames(first_assay) == rownames(rowData))) {
    stop(
      "Names (likely protein names) do not match between the Assay rows and rowData rows"
    )
  }
  new_OlinkAssay(
    names = rownames(assays@data[[1]]),
    assays = assays,
    rowData = rowData,
    colData = colData,
    metadata = metadata,
    negativeControls = negativeControls,
    plateControls = plateControls
  )
}

#' @title OlinkAssayFromNPX
#' @description Convert data from an NPX parquet file to an [`OlinkAssay`] object
#'
#' @param npxData a [`tibble`][`tibble::tibble`] with NPX data
#' @param colData a [`tibble`][`tibble::tibble`] with extra sample metadata to store
#'  in the `colData` slot.
#' @param assay_tables A list of columns in `npxData` to store
#'  in the `assays` slot. (Default: "Count", "ExtNPX", "PCNormalizedNPX", "NPX")
#' @param verbose Should extra information about conversion progress be shown? (Default: FALSE)
#'
#' @export
#' @returns an [`OlinkAssay`][`OlinkAssay::OlinkAssay`] object
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame SimpleList
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

  .metadata[["PlateRef"]] <- stringr::str_split_i(
    string = npxData[["PlateID"]],
    pattern = "_",
    i = 1
  ) |>
    unique()
  # should these, along with the plate and negative controls below,
  # be stored in slots instead? Seems like I'm just dumping a lot in `metadata`
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

  data_cols <- c(
    "Count",
    "ExtNPX",
    "NPX",
    "PCNormalizedNPX",
    "ExtNPX_Corrected",
    "LogProtExp",
    "LogProtExp_Raw"
  )

  # TODO: Should these include the assay controls? E.g. "Extension control 1", "Incubation control 1"...
  negative_control_data <-
    dplyr::filter(
      .data = npxData,
      SampleType == "NEGATIVE_CONTROL"
    ) |>
    dplyr::select(
      Assay,
      OlinkID,
      PlateID,
      tidyselect::any_of(data_cols)
    ) |>
    dplyr::mutate(
      Assay = as.factor(Assay), # not strictly necessary, but given that there's a set number of Assay names, seems appropriate
      PlateID = stringr::str_remove(string = PlateID, pattern = "_plate[0-9]$")
    ) |>
    dplyr::rename(PlateRef = PlateID) |>
    S4Vectors::DataFrame()

  plate_control_data <-
    dplyr::filter(
      .data = npxData,
      SampleType == "PLATE_CONTROL"
    ) |>
    dplyr::select(
      Assay,
      OlinkID,
      PlateID,
      tidyselect::any_of(data_cols)
    ) |>
    dplyr::mutate(
      Assay = as.factor(Assay),
      PlateID = stringr::str_remove(string = PlateID, pattern = "_plate[0-9]$")
    ) |>
    dplyr::rename(PlateRef = PlateID) |>
    S4Vectors::DataFrame()

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
    .coldata <- dplyr::left_join(
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

  OlinkAssay(
    assays = SummarizedExperiment::Assays(sample_df),
    rowData = .rowdata,
    colData = .coldata,
    metadata = .metadata,
    negativeControls = negative_control_data,
    plateControls = plate_control_data
  )
}

#### object validation ####
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

#### assay getter ####
#' @title assay
#' @description Get assay data for OlinkAssay object
#' @details Works the same as that for [`SummarizedExperiment::assay`] except that
#'  if a value for `i` is not provided, [`OlinkAssay::assay`] attempts to extract
#'  the assay named (in decending order) "ExtNPX_Corrected", "ExtNPX",
#'  "PCNormalizedNPX", "Correction", "NPX", "Count"
#'
#' @param i name of the assay in `assays(x)` to retreive.
#' @param withDimnames retain row and column names? (Default: TRUE)
#'
#' @inheritParams SummarizedExperiment::assay
#' @export
setGeneric("assay", function(x, ...) standardGeneric("assay"))

#' @rdname assay
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
          "ExtNPX_Corrected",
          "ExtNPX",
          "LogProtExp",
          "PCNormalizedNPX",
          "NPX",
          "Count"
        ),
        assayNames(x)
      )[[1]]
    } else if (!i %in% assayNames(x)) {
      stop(glue::glue(
        "`{i}` is not a valid assay name this object! Please select one of {glue::glue_collapse(assayNames(oa), sep = ', ')}"
      ))
    }

    SummarizedExperiment::assay(x = x, i = i, withDimnames = withDimnames)
  }
)

#### colData ####
#' @title colData
#'
#' @inheritParams SummarizedExperiment::colData
#' @export
setGeneric("colData", function(x, ...) standardGeneric("colData"))

#' @rdname colData
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

#### rowData ####
#' @title rowData
#'
#' @inheritParams SummarizedExperiment::rowData
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

#### medianCorrection ####
#' @title medianCorrection
#' @description Perform batch correction using per-assay plate assay medians.
#' @details During import of NPX data or when creating an `OlinkAssay`, we calculate
#'  the assay medians for the plate controls (those with their `SampleType` labeled
#'  as "PLATE_CONTROL") and the assay medians and for all samples (`SampleType` ==
#'  "SAMPLE") and store them in the `metadata` slot as "AssayMedians" and "GlobalMedians",
#'  respectively. Here, we use those to calculate a correction factor:
#'  \deqn{
#'  Correction_{assay} = Median_{plate} - \frac{\sum_{}^{plates}median_{plate}}{n_{plates}}
#'  }
#'  \deqn{ExtNPX\_Corrected = ExtNPX - Correction_{assay}}
#'  The corrected values are then stored as "ExtNPX_Corrected" in the assays slot
#'
#' @param method Method to use when correcting values, either "median" or "global median".
#'  Generally one should use the "median" methods; "global median" is more appropriate for
#'  plates where the sample groups were not well randomized among wells.
#' @export
setGeneric("medianCorrection", function(x, ...) {
  standardGeneric("medianCorrection")
})

# TODO: either here or in the concat method we need to check that each plate has
# the necessary AssayMedian/GlobalMedian information
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

    assay(x, i = "ExtNPX_Corrected") <- data_correction |>
      dplyr::select(SampleID, assay, ExtNPX_Corrected) |>
      dplyr::distinct() |>
      tidyr::pivot_wider(
        names_from = "SampleID",
        values_from = "ExtNPX_Corrected"
      ) |>
      tibble::column_to_rownames("assay")

    metadata(x)[["BatchCorrectionMethod"]] <- method

    for (i in c("negativeControls", "plateControls")) {
      slot(x, i) <- merge(
        slot(x, i),
        meds_correction,
        by = c("OlinkID", "PlateRef")
      )
      slot(x, i)[["ExtNPX_Corrected"]] <- slot(x, i)[["ExtNPX"]] -
        slot(x, i)[["Correction"]]
      slot(x, i)[["LogProtExp_Raw"]] <- slot(x, i)[[
        "LogProtExp"
      ]] <- slot(x, i)[["ExtNPX_Corrected"]] + log2(1e5)
    }

    x
  }
)

#### concat ####
#' @title concat
#' @description Concatenate or combine two existing OlinkAssay objects.
#' @details Generally, an OlinkAssay object should initially be created from each
#'  run plate separately so that the factors required for batch correction are
#'  appropriately calculated. Using `concat` to combine the objects ensures that
#'  those factors are retained and, if "median" or "global median" are passed to
#'  `batch_correction_method`, used to correct the ExtNPX values.
#'
#' @param x First object to combine
#' @param y Second object
#' @param batch_correction_method method to use in performing batch correction.
#'  See \link{medianCorrection} for details. Default: NULL
#' @export
setGeneric("concat", function(x, y, ...) standardGeneric("concat"))

#' @export
#' @importFrom S4Vectors mcols
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment Assays
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
      assays = SummarizedExperiment::Assays(merged_assays),
      colData = concat_coldata,
      rowData = `rownames<-`(concat_rowdata, NULL),
      metadata = md_list,
      negativeControls = rbind(x@negativeControls, y@negativeControls),
      plateControls = rbind(x@plateControls, y@plateControls)
    )

    if (!is.null(batch_correction_method)) {
      concat_oa <- medianCorrection(concat_oa, method = batch_correction_method)
    }

    concat_oa
  }
)
