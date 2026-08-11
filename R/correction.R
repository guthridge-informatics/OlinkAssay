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
ctrl_ref <- function(data) {
  data |>
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
global_ref <- function(data) {
  data |>
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
#' @importFrom dplyr left_join select mutate c_across transmute rowwise
#' @importFrom purrr map map2 reduce
#' @importFrom tidyselect contains
#'
#' @export
#' @examples
median_correction <- function(data, meds) {
  # calculate the medians for each 384-well plate
  ref_med <-
    purrr::reduce(
      .x = meds,
      .f = \(x, y) dplyr::left_join(x, y, by = "OlinkID")
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
        dplyr::left_join(x = x, y = ref_med) |>
          dplyr::mutate(Correction = Median - ReferenceMedian)
      }
    )

  data_correction <-
    purrr::map2(
      .x = exdata,
      .y = meds_correction,
      .f = \(x, y) {
        dplyr::left_join(x = x, y = y, by = "OlinkID") |>
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
batch_correction <- function(data, method = c("median", "global median")) {
  method <- match.arg(method)
  if (!is.vector(data)) {
    data <- list(data)
  }
  if (method == "median") {
    meds <- purrr::map(.x = data, .f = ctrl_ref)
  } else if (method == "global median") {
    meds <- purrr::map(.x = data, .f = global_ref)
  }

  median_correction(data, meds)
}
