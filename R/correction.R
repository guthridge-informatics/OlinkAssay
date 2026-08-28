#' @title ctrl_ref
#' @description Calculate the default plate control batch correction
#'
#' @param .data [`tibble::tibble`] with Olink data
#'
#' @returns
#'
#' @importFrom dplyr filter group_by summarize
#' @importFrom stats na.omit
#'
#' @export
#' @examples
ctrl_ref <- function(.data) {
  plate_ref <- stringr::str_split_i(
    string = .data[["PlateID"]],
    pattern = "_",
    i = 1
  ) |>
    unique()

  dplyr::filter(
    .data = .data,
    SampleType == "PLATE_CONTROL" & AssayType == "assay"
  ) |>
    dplyr::group_by(OlinkID) |>
    dplyr::summarise(
      Median = median(stats::na.omit(ExtNPX)),
      Variance = var(stats::na.omit(ExtNPX))
    ) |>
    dplyr::mutate(PlateRef = plate_ref)
}

#' @title global_ref
#' @description To be used for Studies where proper randomization can be verified or where plate controls cannot be used for normalization
#'
#' @param .data [`tibble::tibble`] with Olink data
#'
#' @returns
#'
#' @importFrom dplyr filter group_by summarise
#'
#' @export
#' @examples
global_ref <- function(.data) {
  plate_ref <- stringr::str_split_i(
    string = .data[["PlateID"]],
    pattern = "_",
    i = 1
  ) |>
    unique()

  dplyr::filter(
    .data = .data,
    SampleType == "SAMPLE" & AssayType == "assay"
  ) |> # filter down to just the samples and the assay
    dplyr::group_by(OlinkID) |> # grouped it by just the OlinkID
    dplyr::summarise(
      Median = median(stats::na.omit(ExtNPX)),
      Variance = var(stats::na.omit(ExtNPX))
    ) |>
    dplyr::mutate(PlateRef = plate_ref)
}

#' @title median_correction
#' @description This function takes either global or plate control median frames
#'  and applies differences
#'
#' @param .data [`tibble`][`tibble::tibble`] or list of `tibbles` with Olink data
#' @param medians [`tibble`][`tibble::tibble`] with three columns
#' 1. OlinkID
#' 2. Median
#' 3. Variance
#' 4. PlateRef
#'
#' Generally best to just use the output from `ctrl_ref()` or `global_ref()`.
#'
#' @returns
#'
#' @importFrom dplyr left_join select mutate c_across transmute rowwise join_by
#' @importFrom purrr map map2 reduce
#' @importFrom tidyselect contains
#'
#' @export
#' @examples
median_correction <-
  function(
    .data,
    ...
  ) {
    UseMethod("median_correction")
  }


#' @rdname median_correction
#' @method median_correction tbl_df
#' @exportS3Method oinkqc::median_correction
#' @returns
median_correction.tbl_df <- function(.data, medians) {
  # calculate the medians for each 384-well plate
  ref_med <-
    medians |>
    # calculate the mean of all medians to scale by
    # and make a tibble that contains the OlinkID and the reference median
    dplyr::rowwise() |>
    dplyr::transmute(
      OlinkID = OlinkID,
      ReferenceMedian = mean(dplyr::c_across(tidyselect::contains(
        "median",
        ignore.case = TRUE
      ))),
    )

  meds_correction <-
    dplyr::left_join(x = medians, y = ref_med, by = dplyr::join_by(OlinkID)) |>
    dplyr::mutate(Correction = Median - ReferenceMedian)

  data_correction <-
    dplyr::left_join(
      x = .data,
      y = meds_correction,
      by = dplyr::join_by(OlinkID)
    ) |>
    dplyr::mutate(ExtNPX_Corrected = ExtNPX - Correction)

  data_correction
}

#' @rdname median_correction
#' @method median_correction list
#' @export
#' @returns
median_correction.list <- function(.data, medians) {
  # calculate the medians for each 384-well plate
  ref_med <-
    purrr::reduce(
      .x = medians,
      .f = \(x, y) dplyr::left_join(x, y, by = dplyr::join_by(OlinkID))
    ) |> # left join all the run median and variance per each run
    # calculate the mean of all medians to scale by
    # and make a tibble that contains the OlinkID and the reference median
    dplyr::rowwise() |>
    dplyr::transmute(
      OlinkID = OlinkID,
      ReferenceMedian = mean(dplyr::c_across(tidyselect::contains("Median")))
    )

  meds_correction <-
    purrr::map(
      .x = medians,
      .f = \(x) {
        dplyr::left_join(x = x, y = ref_med, by = dplyr::join_by(OlinkID)) |>
          dplyr::mutate(Correction = Median - ReferenceMedian)
      }
    )

  data_correction <-
    purrr::map2(
      .x = .data,
      .y = meds_correction,
      .f = \(x, y) {
        dplyr::left_join(x = x, y = y, by = dplyr::join_by(OlinkID)) |>
          dplyr::mutate(ExtNPX_Corrected = ExtNPX - Correction)
      }
    )

  data_correction
}

#' @title batch_correction
#' @description Main routing function for the batch correction and normalization
#'
#' @param .data [`tibble::tibble`] with Olink data
#' @param method method to use with performing batch correction one of "median"
#'  or "global median". (Default: "median")
#'
#' @returns
#'
#' @importFrom purrr map
#'
#' @export
#' @examples
batch_correction <-
  function(
    .data,
    ...
  ) {
    UseMethod("batch_correction")
  }

#' @rdname batch_correction
#' @method batch_correction tbl_df
#' @exportS3Method olinkqc::batch_correction
#' @returns
batch_correction.tbl_df <- function(
  .data,
  method = c("median", "global median")
) {
  method <- match.arg(method)
  # if (!is.vector(.data)) {
  #   .data <- list(.data)
  # }
  if (method == "median") {
    meds <- ctrl_ref(.data)
  } else if (method == "global median") {
    meds <- global_ref(.data)
  }

  median_correction(.data, meds)
}

#' @rdname batch_correction
#' @method batch_correction list
#' @export
#' @returns
batch_correction.list <- function(
  .data,
  method = c("median", "global median")
) {
  method <- match.arg(method)

  if (method == "median") {
    meds <- purrr::map(.x = .data, .f = ctrl_ref)
  } else if (method == "global median") {
    meds <- purrr::map(.x = .data, .f = global_ref)
  }

  median_correction(.data, meds)
}
