#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.OlinkAssay <- setClass("OlinkAssay", contains = "SummarizedExperiment")

#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
OlinkAssay <- function(.data, ...) {
  se <- as_SummarizedExperiment(.data, ...)
  .OlinkAssay(se)
}

#' @title assay
#' @description Get assay data for OlinkAssay object
#'
#' @export
setGeneric("assay", function(x, ...) standardGeneric("assay"))

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
          "PCNormalizedNPX",
          "Correction",
          "NPX",
          "Count"
        ),
        colnames(rowData(x))
      )[[1]]
    }
    SummarizedExperiment::assay(x = x, i = i, withDimnames = withDimnames)
  }
)

#' @export
setGeneric("npx", function(x, ...) standardGeneric("npx"))

#' @export
#' @importFrom SummarizedExperiment assay
setMethod(
  f = "npx",
  signature = "OlinkAssay",
  definition = function(x, withDimnames = TRUE) {
    SummarizedExperiment::assay(x = x, i = "NPX", withDimnames = withDimnames)
  }
)


#' @export
setGeneric("colData", function(x, ...) standardGeneric("colData"))

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
