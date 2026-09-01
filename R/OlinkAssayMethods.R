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
      tibble::rownames_to_column("Assay") |>
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
    colnames(sample_assay_qc) <- colnames(sample_assay_qc) |>
      stringr::str_remove(pattern = "^X")
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
    colnames(new_logprotext) <- colnames(new_logprotext) |>
      stringr::str_remove(pattern = "^X")
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
