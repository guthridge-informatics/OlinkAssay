#### scaleNPX ####
#' @title scaleNPX
#' @description Scale ExtNPX data by adding \eqn{\text{log}_2(10^{5})} to create a non-negative
#' data matrix
#'
#' @param x
#' @export
setGeneric("scaleNPX", function(x, ...) standardGeneric("scaleNPX"))

#' @export
#' @importFrom SummarizedExperiment assay assay<- assayNames
setMethod(
  f = "scaleNPX",
  signature = "OlinkAssay",
  definition = function(x) {
    if ("ExtNPX_Corrected" %in% assayNames(x)) {
      source_assay <- "ExtNPX_Corrected"
    } else {
      warning("Batch corrected data not found - using `ExtNPX` assay instead.")
      source_assay <- "ExtNPX"
    }
    assay(x, "LogProtExp") <- assay(x, source_assay) + log2(1e5)
    assay(x, "LogProtExp_Raw") <- assay(x, source_assay) + log2(1e5)
    x
  }
)

#### level2Process ####
#' @title level2Process
#' @description Calculate assay and sample QC metrics matching those as described by
#'  the Olink HT data standard document's level 2 requirements
#'
#' @param x
#' @export
setGeneric("level2Process", function(x, ...) standardGeneric("level2Process"))

#' @export
#' @importFrom methods slot
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom S4Vectors DataFrame
#' @importFrom tibble as_tibble rownames_to_column column_to_rownames
#' @importFrom dplyr group_by summarize mutate if_else select left_join join_by case_when distinct n
#' @importFrom stats median na.omit quantile sd
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tidyselect where
#' @importFrom stringr str_remove
#'
setMethod(
  f = "level2Process",
  signature = "OlinkAssay",
  definition = function(x) {
    if (!"LogProtExp" %in% assayNames(x)) {
      stop(
        "`LogProtExp` was not found in the assay slot. Please run `scaleNPX` first."
      )
    }
    ht_nc_vals <-
      slot(x, "negativeControls") |>
      tibble::as_tibble() |>
      dplyr::group_by(Assay, OlinkID) %>%
      dplyr::summarise(
        median_nc = stats::median(stats::na.omit(LogProtExp)),
        iqr_nc = as.numeric(stats::quantile(stats::na.omit(LogProtExp), 0.75)),
        .groups = 'drop'
      )

    ht_pc_vals <-
      slot(x, "plateControls") |>
      tibble::as_tibble() |>
      dplyr::group_by(Assay, OlinkID) %>%
      dplyr::summarise(
        pc_cv = 100 * stats::sd(LogProtExp) / mean(LogProtExp),
        .groups = 'drop'
      ) %>%
      dplyr::mutate(
        high_var_assay = dplyr::if_else(
          condition = pc_cv >= 20,
          true = "High Variance",
          false = "Pass"
        )
      ) |>
      dplyr::select(-pc_cv)

    ht_scaled_npx_sample <-
      SummarizedExperiment::assay(x, "LogProtExp") |>
      tibble::as_tibble(rownames = "Assay") |>
      tidyr::pivot_longer(
        cols = tidyselect::where(is.numeric),
        names_to = "SampleID",
        values_to = "LogProtExp"
      ) |>
      dplyr::left_join(ht_nc_vals, by = dplyr::join_by("Assay")) |>
      dplyr::left_join(ht_pc_vals, by = dplyr::join_by("Assay", "OlinkID")) |>
      dplyr::mutate(
        sample_level_qc = dplyr::case_when(
          LogProtExp < median_nc ~ "Below LLOD",
          LogProtExp < iqr_nc ~ "Below LLOQ",
          TRUE ~ "Pass"
        ),
        LogProtExp = dplyr::case_when(
          LogProtExp < median_nc ~ 0,
          LogProtExp < iqr_nc ~ iqr_nc,
          TRUE ~ LogProtExp
        )
      )

    sample_assay_qc <- ht_scaled_npx_sample |>
      dplyr::select(SampleID, Assay, sample_level_qc) |>
      dplyr::distinct() |>
      tidyr::pivot_wider(
        names_from = "SampleID",
        values_from = "sample_level_qc"
      ) |>
      tibble::column_to_rownames("Assay") |>
      S4Vectors::DataFrame()
    SummarizedExperiment::assay(x, "SampleAssayQC") <- sample_assay_qc

    new_logprotext <-
      ht_scaled_npx_sample |>
      dplyr::select(SampleID, Assay, LogProtExp) |>
      tidyr::pivot_wider(
        names_from = "SampleID",
        values_from = "LogProtExp"
      ) |>
      tibble::column_to_rownames("Assay") |>
      S4Vectors::DataFrame()
    SummarizedExperiment::assay(x, "LogProtExp") <- new_logprotext

    ht_scaled_npx_assay <-
      ht_scaled_npx_sample |>
      dplyr::select(SampleID, Assay, sample_level_qc) |>
      dplyr::distinct() |>
      dplyr::group_by(Assay, sample_level_qc) |>
      dplyr::summarise(
        percentage = 100 * dplyr::n() / nrow(ht_scaled_npx_sample),
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
        ),
        `Below LLOQ` = `Below LLOQ` + `Below LLOD`,
        assay_level_qc = dplyr::case_when(
          `Below LLOD` > 75 ~ "Categorical",
          `Below LLOQ` > 50 ~ "Semi-Continuous",
          T ~ "Continuous"
        )
      ) |>
      dplyr::select(Assay, Pass, assay_level_qc)

    new_rowdata <-
      ht_scaled_npx_sample |>
      dplyr::select(Assay, median_nc, iqr_nc, high_var_assay) |>
      dplyr::distinct() |>
      dplyr::left_join(
        tibble::as_tibble(rowData(x), rownames = "Assay"),
        by = dplyr::join_by("Assay")
      ) |>
      dplyr::left_join(
        ht_scaled_npx_assay,
        by = dplyr::join_by("Assay")
      ) |>
      tibble::column_to_rownames("Assay") |>
      S4Vectors::DataFrame()

    rowData(x) <- new_rowdata

    x
  }
)

#### batchCorrection ####
# TODO: either here or in the concat method we need to check that each plate has
# the necessary AssayMedian/GlobalMedian information
#' @importFrom tibble as_tibble column_to_rownames
#' @importFrom dplyr select rowwise transmute c_across left_join join_by mutate distinct
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom tidyselect contains where
#' @importFrom S4Vectors metadata merge
#' @importFrom stringr str_remove
#' @importFrom SummarizedExperiment assay
#' @importFrom methods slot
#' @export
#'
setMethod(
  f = "batchCorrection",
  signature = "OlinkAssay",
  definition = function(
    x,
    method = c("median", "global median")
  ) {
    correction_values <- switch(
      EXPR = method,
      "median" = S4Vectors::metadata(x)[["AssayMedians"]],
      "global median" = S4Vectors::metadata(x)[["GlobalMedians"]]
    )

    median_means <-
      correction_values |>
      tibble::as_tibble() |>
      dplyr::select(-Variance) |>
      tidyr::pivot_wider(
        names_from = "PlateRef",
        names_glue = "{PlateRef}_median",
        # names_from = "PlateID",
        # names_glue = "{PlateID}_median",
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
    # dplyr::select(SampleID, PlateID)

    npx_values <- SummarizedExperiment::assay(x, i = "ExtNPX") |>
      tibble::as_tibble(rownames = "assay") |>
      tidyr::pivot_longer(
        tidyselect::where(is.numeric),
        names_to = "SampleID",
        values_to = "ExtNPX"
      )

    corrected_data <-
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
        # by = dplyr::join_by(OlinkID, PlateID)
      ) |>
      dplyr::mutate(ExtNPX_Corrected = ExtNPX - Correction)

    SummarizedExperiment::assay(x, i = "ExtNPX_Corrected") <- corrected_data |>
      dplyr::select(SampleID, assay, ExtNPX_Corrected) |>
      dplyr::distinct() |>
      tidyr::pivot_wider(
        names_from = "SampleID",
        values_from = "ExtNPX_Corrected"
      ) |>
      tibble::column_to_rownames("assay") |>
      as.matrix()

    S4Vectors::metadata(x)[["BatchCorrectionMethod"]] <- method

    for (i in c("negativeControls", "plateControls")) {
      methods::slot(x, i) <- S4Vectors::merge(
        methods::slot(x, i),
        meds_correction,
        by = c("OlinkID", "PlateRef")
        # by = c("OlinkID", "PlateID")
      )
      methods::slot(x, i)[["ExtNPX_Corrected"]] <- methods::slot(x, i)[[
        "ExtNPX"
      ]] -
        methods::slot(x, i)[["Correction"]]
      methods::slot(x, i)[["LogProtExp_Raw"]] <- methods::slot(x, i)[[
        "LogProtExp"
      ]] <- methods::slot(x, i)[["ExtNPX_Corrected"]] + log2(1e5)
    }
    x
  }
)

#' @title SampleAssayQC
#' @description Calculate Assay and Sample QC stats
#' @param x OlinkAssay object
#' @export
setGeneric("SampleAssayQC", function(x, ...) {
  standardGeneric("SampleAssayQC")
})

#' @export
setMethod(
  f = "SampleAssayQC",
  signature = "OlinkAssay",
  definition = function(x) {
    assay_tbl <- rowData(x) |>
      tibble::as_tibble(rownames = "Assay") |>
      tidyr::unnest(PerPlateAssayQC) |>
      tidyr::separate_wider_delim(
        PerPlateAssayQC,
        delim = ":",
        names = c("PlateID", "AssayQC")
      ) |>
      dplyr::group_by(PlateID, OlinkID, AssayQC) |> # group by PlateID and OlinkID to get the split among the Assay QC
      dplyr::tally() |> # tabulate the counts
      dplyr::ungroup() |> # this is to allow the ungrouping of the tibble and return to orignal tibble without grouping
      dplyr::group_by(PlateID, AssayQC) |> # grouping to allow tabulation of assays
      dplyr::tally() |>
      dplyr::ungroup()

    sample_tbl <-
      colData(x) |>
      tibble::as_tibble(rownames = "SampleID") |>
      tidyr::unnest(PerAssaySampleQC) |>
      tidyr::separate_wider_delim(
        PerAssaySampleQC,
        delim = ":",
        names = c("Assay", "SampleQC")
      ) |>
      dplyr::group_by(PlateID, SampleID, SampleQC) |> # group by PlateID and OlinkID to get the split among the Assay QC
      dplyr::tally() |> # tabulate the counts
      dplyr::mutate(Frequency = prop.table(n)) |> # calculate the frequencies within each sample
      dplyr::ungroup()

    list(
      AssayOlinkQC = assay_tbl,
      SampleOlinkQC = sample_tbl
    )
  }
)

#' @title retrieveSampleQC
#' @description Extract the SampleQC results from colData and reformat
#' @param x OlinkAssay object
#' @export
setGeneric("retrieveSampleQC", function(x, ...) {
  standardGeneric("retrieveSampleQC")
})

#' @export
setMethod(
  f = "retrieveSampleQC",
  signature = "OlinkAssay",
  definition = function(x) {
    colData(x) |>
      tibble::as_tibble(rownames = "SampleID") |>
      tidyr::unnest(PerAssaySampleQC) |>
      tidyr::separate_wider_delim(
        cols = PerAssaySampleQC,
        delim = ":",
        names = c("Assay", "SampleQC")
      ) |>
      dplyr::select(SampleID, PlateID, Assay, SampleQC)
  }
)

#' @title retrieveAssayQC
#' @description Extract the AssayQC results from rowData and reformat
#' @param x OlinkAssay object
#' @export
setGeneric("retrieveAssayQC", function(x, ...) {
  standardGeneric("retrieveAssayQC")
})

#' @export
setMethod(
  f = "retrieveAssayQC",
  signature = "OlinkAssay",
  definition = function(x) {
    rowData(x) |>
      tibble::as_tibble(rownames = "Assay") |>
      tidyr::unnest(c(PerPlateAssayQC, PerPlateAssayWarn)) |>
      tidyr::separate_wider_delim(
        PerPlateAssayQC,
        delim = ":",
        names = c("PlateID", "AssayQC")
      ) |>
      tidyr::separate_wider_delim(
        PerPlateAssayWarn,
        delim = ":",
        names = c(NA, "AssayQCWarn")
      )
  }
)
