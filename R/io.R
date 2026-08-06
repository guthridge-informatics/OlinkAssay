#' @title multifile_read
#' @description Read multiple files into a list and name them based on the file
#' name minus the file extension
#'
#' @param directory path to where the individual Olink run data can be found
#' @param file_extension what kind of file should be loaded: `"parquet"` or `"excel"` (default `"parquet"`)
#'
#' @importFrom stringr str_split_i str_remove
#' @importFrom readxl read_excel
#' @importFrom readr read_csv
#' @importFrom purrr map
#' @importFrom arrow read_parquet
#'
#' @returns A named list of tibbles
#'
.reader_func <- function(directory, file_extension = c("parquet", "excel")) {
  file_extension <- match.arg(file_extension, several.ok = FALSE)
  # filter on the regex ^[[:alnum:]] to stop from loading
  # Excel's temporary files

  if (file_extension == "parquet") {
    pattern <- "\\.parquet$"
    read_func <- arrow::read_parquet
  } else if (file_extension == "excel") {
    pattern <- "^[[:alnum:]].*\\.xlsx"
    read_func <- \(x) {
      {
        readxl::read_excel(
          path = x,
          sheet = 1
        )
      } |>
        dplyr::select(sample_id = Assay_Sample_ID, Project)
    }
  } else if (file_extension == "csv") {
    pattern <- "^[[:alnum:]].*\\.csv"
    read_func <- \(x) {
      readr::read_csv(
        file = x,
        col_select = c(sample_id = Assay_Sample_ID, Project)
      )
    }
  } else {
    stop("File must be either in parquet, excel, or comma-delimited format")
  }

  files_run <- list.files(
    path = directory,
    pattern = pattern,
    full.names = TRUE
  ) |>
    na.omit() |>
    as.character()
  if (length(files_run) < 1) {
    stop(stringr::str_glue("No {file_extension} files were found!"))
  }
  data <- purrr::map(.x = files_run, .f = read_func)

  names(data) <- stringr::str_split_i(
    string = files_run,
    pattern = .Platform$file.sep,
    i = -1
  ) |>
    stringr::str_remove(pattern = paste(".", file_extension, sep = ""))

  data
}

#' @title ingest_olink_data
#' @description Read multiple Olink output files into a list and name them based on the file name
#'
#' @param directory path to where the individual Olink run data can be found
#'
#' @returns A named list of tibbles
#' @export
#'
ingest_olink_data <- function(directory) {
  .reader_func(directory, file_extension = "parquet")
}

#' @title ingest_manifest
#' @description Read in the manifest of samples for an Olink run \
#'
#' @param directory path to where Olink run manifest
#'
#' @returns A named list of tibbles
#' @export
#'
ingest_manifest <- function(directory) {
  .reader_func(directory, file_extension = "excel")
}


#' @title multifile_write
#' @description Write multiple files into a list and name them based on the file name without file extension
#'
#' @param data A named list of [tibble::tibble]
#' @param file_extension Output type: `"parquet"` or `"csv"` (default `"parquet"`)
#'
#' @importFrom purrr imap
#' @importFrom arrow write_parquet
#' @importFrom readr write_csv
#'
#' @returns Nothing
#'
#' @export
#' @examples
multifile_write <- function(data, file_extension = c("parquet", "csv")) {
  file_extension <- match.arg(file_extension)

  if (file_extension == "parquet") {
    purrr::imap(
      .x = data,
      .f = \(x, name) {
        arrow::write_parquet(
          x,
          sink = paste0("SDP/Level_1/", name, ".parquet")
        )
      }
    )
  } else if (file_extension == "csv") {
    purrr::imap(
      .x = data,
      .f = \(x) {
        readr::write_csv(
          x = x,
          file = paste0("SDP/Level_1/", x, ".csv")
        )
      }
    )
  }
}


#' @title Olink_Reader
#' @description Setup `"directory"` according to the project format described in the Olink data standards
#' document, and then read both parquet files and manfests
#' and reads them into a list format
#'
#' @param directory project directory to which `"code"`, `"level 1"`, and `"level 2"` files should be
#' written and where the output from NPXExplorer and the sample manifests can be found
#'
#' @importFrom purrr map
#' @returns list with `"data"` containing a list of [`tibbles`][`tibble::tibble`]
#'
#' @export
#' @examples
Olink_Reader <- function(directory) {
  # Defining the directory all Olink data is stored in

  # make sub-directory for SDP hiearchical file structure
  purrr::map(
    .x = c(
      root = paste0(directory, "/SDP"), # root directory for SDP
      lvl1 = paste0(directory, "/SDP/Level_1"), # level 1 directory for SDP
      lvl2 = paste0(directory, "/SDP/Level_2"), # level 2 directory for SDP
      code = paste0(directory, "/SDP/code")
    ),
    .f = dir.create,
    recursive = TRUE
  )

  # Return parquet files
  list(
    data = ingest_olink_data(directory),
    manifest = ingest_manifest(directory)
  )
}
