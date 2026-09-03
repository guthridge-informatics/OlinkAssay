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
#' @param extra_metadata A named list of items to add to the metadata list
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
  extra_metadata = NULL,
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
    message("Extracting plate ref")
  }
  plate_ref <- stringr::str_remove(
    string = npxData[["PlateID"]],
    pattern = "_.*$"
  ) |>
    unique()

  if (verbose) {
    message("Sanitizing sample IDs")
  }
  # we have to sanitize the sampleID here
  # because if we don't, S4Vectors::DataFrame *will*
  # which will result in mismatched sampleIDs

  # ideally, I'd use janitor::make_clean_names, but
  # my gods it was slow
  npxData[["SampleID"]] <- npxData[["SampleID"]] |>
    stringr::str_replace("^([0-9])", "X\\1") |>
    stringr::str_replace("[\\|\\.-]", "_")

  if (verbose) {
    message("Extracting metadata")
  }
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

  if (!is.null(extra_metadata)) {
    .metadata <- c(.metadata, extra_metadata)
  }

  data_cols <- c(
    "Count",
    "ExtNPX",
    "NPX",
    "PCNormalizedNPX",
    "ExtNPX_Corrected",
    "LogProtExp",
    "LogProtExp_Raw",
    "AssayQC",
    "SampleQC"
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
          tibble::column_to_rownames(var = "Assay") |>
          as.matrix()
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
      Assay,
      SampleQC,
      OSICategory,
      OSISummary,
      OSITimeToCentrifugation,
      OSIPreparationTemperature
    ) |>
    dplyr::distinct() |>
    dplyr::mutate(
      PerAssaySampleQC = list(rlang::set_names(
        glue::glue("{Assay}:{SampleQC}"),
        Assay
      )),
      .by = SampleID
    ) |>
    dplyr::select(-Assay, -SampleQC) |>
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
      Normalization,
      Block,
      PlateID,
      AssayQC,
      AssayQCWarn
    ) |>
    dplyr::distinct() |>
    dplyr::mutate(
      PerPlateAssayQC = list(rlang::set_names(
        glue::glue("{PlateID}:{AssayQC}"),
        PlateID
      )),
      PerPlateAssayWarn = list(rlang::set_names(
        glue::glue("{PlateID}:{AssayQCWarn}"),
        PlateID
      )),
      .by = Assay
    ) |>
    dplyr::select(
      -PlateID,
      -AssayQC,
      -AssayQCWarn
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

#' @title OlinkAssayFromDisk
#' @export
#' @importFrom arrow read_parquet
OlinkAssayFromDisk <- function(
  npx_file = NULL,
  metadata_file = NULL,
  metadata_sheet = 1,
  sample_column = "SampleID",
  project_column = "Project",
  additional_columns = NULL,
  extra_metadata = NULL,
  verbose = FALSE
) {
  if (is.character(npx_file)) {
    if (verbose) {
      message("Loading NPX file")
    }
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

    if (verbose) {
      message("Loading manifest")
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

  se <- OlinkAssayFromNPX(
    npxData = df,
    colData = md,
    extra_metadata = extra_metadata,
    verbose = verbose
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

# Lifted from https://github.com/Bioconductor/SummarizedExperiment/blob/ffe9db3f6666e215296d64dc74e9aae7600bddf0/R/SummarizedExperiment-class.R#L139-L150
setGeneric("colData<-", function(x, ..., value) standardGeneric("colData<-"))

setReplaceMethod(
  "colData",
  c("OlinkAssay", "DataFrame"),
  function(x, ..., value) {
    if (nrow(value) != ncol(x)) {
      stop("nrow of supplied 'colData' must equal ncol of object")
    }
    x <- updateObject(x, check = FALSE)
    # this seems like a not-great idea, but it is what is in the {SummarizedExperiment} package?
    BiocGenerics:::replaceSlots(x, colData = value, check = FALSE)
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

# Taken from https://github.com/Bioconductor/SummarizedExperiment/blob/ffe9db3f6666e215296d64dc74e9aae7600bddf0/R/SummarizedExperiment-class.R#L125
setGeneric("rowData<-", function(x, ..., value) standardGeneric("rowData<-"))

setReplaceMethod("rowData", "OlinkAssay", function(x, ..., value) {
  S4Vectors::`mcols<-`(x, ..., value = value)
})

#### batchCorrection ####
#' @title batchCorrection
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
setGeneric("batchCorrection", function(x, ...) {
  standardGeneric("batchCorrection")
})

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
#'  See \link{batchCorrection} for details. Default: NULL
#' @export
setGeneric("concat", function(x, y, ...) standardGeneric("concat"))

#' @export
#' @importFrom S4Vectors mcols metadata
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
    concat_rowdata <- .combineRowData(rowData(x), rowData(y))

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
        c(S4Vectors::metadata(x)[[z]], S4Vectors::metadata(y)[[z]])
      }
    ) |>
      stats::setNames(single_value_cols)

    md_list[["AssayMedians"]] <- rbind(
      S4Vectors::metadata(x)[["AssayMedians"]],
      S4Vectors::metadata(y)[["AssayMedians"]]
    )
    md_list[["GlobalMedians"]] <- rbind(
      S4Vectors::metadata(x)[["GlobalMedians"]],
      S4Vectors::metadata(y)[["GlobalMedians"]]
    )

    # `base_assays` being those that are present when naive data
    # is first imported. If we are merging data, I am not yet sure
    # that merging already corrected data should be included and am
    # instead relying on repeating any previous corrections.
    # will have to test
    base_assays <- c("Count", "ExtNPX", "PCNormalizedNPX", "NPX")
    # would love to merge the tables here without dplyr, but using
    # `merge` from base/S4Vectors results in a loss of row order;
    # trying to preserve it involves a bit of messier code
    merged_assays <- purrr::map(
      .x = base_assays,
      .f = \(z) {
        dplyr::left_join(
          tibble::as_tibble(
            SummarizedExperiment::assay(x, i = z),
            rownames = "Assay"
          ),
          tibble::as_tibble(
            SummarizedExperiment::assay(y, i = z),
            rownames = "Assay"
          ),
          by = dplyr::join_by("Assay")
        ) |>
          tibble::column_to_rownames("Assay") |>
          as.matrix()
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
      concat_oa <- batchCorrection(concat_oa, method = batch_correction_method)
    }

    concat_oa
  }
)

.combineRowData <- function(t1, t2) {
  purrr::map(
    .x = list(t1, t2),
    .f = \(x) {
      x |>
        tibble::as_tibble(rownames = "Assay") |>
        tidyr::unnest(
          c(
            PerPlateAssayQC,
            PerPlateAssayWarn
          )
        )
    }
  ) |>
    purrr::list_rbind() |>
    tidyr::separate_wider_delim(
      cols = c(PerPlateAssayQC, PerPlateAssayWarn),
      delim = ":",
      names_sep = "_"
    ) |>
    dplyr::mutate(
      PerPlateAssayQC = list(rlang::set_names(
        glue::glue("{PerPlateAssayQC_1}:{PerPlateAssayQC_2}"),
        PerPlateAssayQC_1
      )),
      PerPlateAssayWarn = list(rlang::set_names(
        glue::glue("{PerPlateAssayWarn_1}:{PerPlateAssayWarn_2}"),
        PerPlateAssayWarn_1
      )),
      .by = Assay
    ) |>
    dplyr::select(
      -tidyselect::matches("[0-9]+$")
    ) |>
    dplyr::distinct() |>
    tibble::column_to_rownames(var = "Assay") |>
    S4Vectors::DataFrame()
}
