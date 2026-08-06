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
#' @importFrom dplyr left_join select contains mutate
#' @importFrom magrittr `%>%`
#' @importFrom purrr map map2
#'
#' @export
#' @examples
median_correction <- function(data, meds) {
  # calculate the medians for each 384-well plate
  ref_med <-
    Reduce(
      f = \(x, y) dplyr::left_join(x, y, by = "OlinkID"),
      x = meds
    ) |> # left join all the run median and variance per each run
    dplyr::select(dplyr::contains("Median")) |> # filter only the median column names
    apply(1, mean) %>% # calculate the mean of all medians to scale it to it
    data.frame(OlinkID = meds[[1]]$OlinkID, ReferenceMedian = .) # make a resulting data frame that contains the OlinkID and the reference median

  meds_correction <-
    purrr::map(
      .x = meds,
      .f = \(x) dplyr::left_join(x, ref_med)
    ) %>%
    purrr::map(
      .f = \(x) {
        dplyr::mutate(x, Correction = Median - ReferenceMedian)
      }
    )

  data_correction <-
    purrr::map2(
      .x = data,
      .y = meds_correction,
      .f = \(x, y) dplyr::left_join(x, y, by = "OlinkID")
    ) %>%
    purrr::map(.f = \(x) {
      dplyr::mutate(x, ExtNPX_Corrected = ExtNPX - Correction)
    })

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

  if (method == "median") {
    meds <- purrr::map(.x = data, .F = ctrl_ref)
  } else if (method == "global median") {
    meds <- purrr::map(.x = data, .F = global_ref)
  }

  median_correction(data, meds)
}

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
