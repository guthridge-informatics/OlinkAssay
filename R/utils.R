#' @title .unique_rowvals
#' @description Return the unique number of values in each column
#'
#' @param .data [`tibble::tibble`]
#'
#' @returns [`tibble::tibble`] with one value per column of `.data` with the number
#'  distinct values in each
#'
#' @importFrom dplyr transmute across distinct n_distinct
#' @importFrom tidyselect everything
#' @importFrom tidyr pivot_longer
#' @examples
#' \dontrun{
#' data(mtcars)
#' .unique_rowvals(mtcars)
#' }
.unique_rowvals <- function(.data) {
  .data |>
    dplyr::transmute(
      dplyr::across(
        .cols = tidyselect::everything(),
        .fns = dplyr::n_distinct
      )
    ) |>
    dplyr::distinct() |>
    tidyr::pivot_longer(
      cols = tidyselect::everything(),
      names_to = "column",
      values_to = "n"
    )
}

#' @title .extract_dfs
#' @description Remove controls for either assays, samples, or both from
#'  the initial NPX dataframe
#'
#' @param df
#' @param filter_assay_controls
#' @param filter_sample_controls
#'
#' @returns
#'
#' @examples
.extract_dfs <- function(
  df,
  filter_assay_controls = FALSE,
  filter_sample_controls = FALSE
) {
  dfs <- df |>
    dplyr::filter(
      stringr::str_detect(
        string = AssayType,
        pattern = "ctrl",
        negate = filter_assay_controls
      ) &
        stringr::str_detect(
          string = SampleType,
          pattern = "CONTROL",
          negate = filter_sample_controls
        )
    ) |>
    dplyr::group_by(SampleType) |>
    dplyr::group_split(.keep = FALSE)

  names(dfs) <- df |>
    dplyr::filter(
      stringr::str_detect(
        string = AssayType,
        pattern = "ctrl",
        negate = filter_assay_controls
      ) &
        stringr::str_detect(
          string = SampleType,
          pattern = "CONTROL",
          negate = filter_sample_controls
        )
    ) |>
    dplyr::group_by(SampleType) |>
    dplyr::group_keys() |>
    unlist()
  dfs
}
