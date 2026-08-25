#' @title ctrl_ref
#' @description Calculate the default plate control batch correction
#'
#' @param data
#'
#' @returns
#'
#' @importFrom dplyr filter group_by summarize
#'
#' @export
#' @examples
ctrl_ref <- function(.data) {
  .data |>
    dplyr::filter(SampleType == "PLATE_CONTROL" & AssayType == "assay") |>
    dplyr::group_by(OlinkID) |>
    dplyr::summarise(
      Median = median(na.omit(ExtNPX)),
      Variance = var(na.omit(ExtNPX))
    )
}

#' @title global_ref
#' @description To be used for Studies where proper randomization can be verified or where plate controls cannot be used for normalization
#'
#' @param data
#'
#' @returns
#'
#' @importFrom dplyr filter group_by summarise
#'
#' @export
#' @examples
global_ref <- function(.data) {
  .data |>
    dplyr::filter(SampleType == "SAMPLE" & AssayType == "assay") |> # filter down to just the samples and the assay
    dplyr::group_by(OlinkID) |> # grouped it by just the OlinkID
    dplyr::summarise(
      Median = median(na.omit(ExtNPX)), # calculate the median
      Variance = var(na.omit(ExtNPX))
    ) # calculate variance
}

#' @title median_correction
#' @description This function takes either global or plate control median frames and applies differences
#'
#' @param data
#' @param meds
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
    object,
    ...
  ) {
    UseMethod("median_correction")
  }


#' @rdname median_correction
#' @method median_correction tibble
#' @exportS3Method oinkqc::median_correction
#' @returns
median_correction.tbl_df <- function(.data, meds) {
  # calculate the medians for each 384-well plate
  ref_med <-
    meds |>
    # calculate the mean of all medians to scale by
    # and make a tibble that contains the OlinkID and the reference median
    dplyr::rowwise() |>
    dplyr::transmute(
      OlinkID = OlinkID,
      ReferenceMedian = mean(dplyr::c_across(tidyselect::contains("Median")))
    )

  meds_correction <-
    dplyr::left_join(x = meds, y = ref_med, by = dplyr::join_by(OlinkID)) |>
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
median_correction.list <- function(.data, meds) {
  # calculate the medians for each 384-well plate
  ref_med <-
    purrr::reduce(
      .x = meds,
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
      .x = meds,
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
#' @param data
#' @param method
#'
#' @returns
#'
#' @importFrom purrr map
#'
#' @export
#' @examples
batch_correction <-
  function(
    object,
    ...
  ) {
    UseMethod("batch_correction")
  }

#' @rdname batch_correction
#' @method batch_correction tibble
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
