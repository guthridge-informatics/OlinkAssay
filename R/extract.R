extract_dfs <- function(
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
