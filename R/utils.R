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
