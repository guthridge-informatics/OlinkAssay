#' @title multifile_read
#' @description Read multiple files into a list and name them based on the file
#' name minus the file extension
#'
#' @param input path to where the individual Olink run data can be found
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
.reader_func <- function(
  input,
  file_extension = c("parquet", "excel"),
  sample_column = NULL,
  project_column = NULL,
  manifest_sheet = 1
) {
  if (!is.null(sample_column)) {
    sample_column <- rlang::sym(sample_column)
  }
  if (!is.null(project_column)) {
    project_column <- rlang::sym(project_column)
  }

  file_extension <- match.arg(file_extension, several.ok = FALSE)
  # filter on the regex ^[[:alnum:]] to stop from loading
  # Excel's temporary files

  if (file_extension == "parquet") {
    pattern <- "\\.parquet$"
    ext <- "parquet"
    read_func <- arrow::read_parquet
  } else if (file_extension == "excel") {
    pattern <- "^[[:alnum:]].*\\.xlsx"
    ext <- "xlsx"
    read_func <- \(x) {
      {
        readxl::read_excel(
          path = x,
          sheet = manifest_sheet
        )
      } |>
        dplyr::select(
          SampleID = {{ sample_column }},
          Project = {{ project_column }}
        )
    }
  } else if (file_extension == "csv") {
    pattern <- "^[[:alnum:]].*\\.csv"
    ext <- "csv"
    read_func <- \(x) {
      readr::read_csv(
        file = x,
      ) |>
        dplyr::select(
          SampleID = {{ sample_column }},
          Project = {{ project_column }}
        )
    }
  } else {
    stop("File must be either in parquet, excel, or comma-delimited format")
  }

  files_run <- list.files(
    path = input,
    pattern = pattern,
    full.names = TRUE
  ) |>
    na.omit() |>
    as.character()
  if (length(files_run) < 1) {
    stop(stringr::str_glue("No {file_extension} files were found!"))
  }
  .data <- purrr::map(.x = files_run, .f = read_func)

  names(.data) <- stringr::str_split_i(
    string = files_run,
    pattern = .Platform$file.sep,
    i = -1
  ) |>
    stringr::str_remove(pattern = paste(".", ext, sep = ""))

  .data
}

#' @title ingest_olink_data
#' @description Read multiple Olink output files into a list and name them based on the file name
#'
#' @param input path to where the individual Olink run data can be found
#'
#' @returns A named list of tibbles
#' @export
#'
ingest_olink_data <- function(input) {
  .reader_func(input, file_extension = "parquet")
}

#' @title ingest_manifest
#' @description Read in the manifest of samples for an Olink run \
#'
#' @param input path to where Olink run manifest. Manifest MUST BE in excel format
#' and MUST HAVE a `SampleID` and `project` columns.
#' @param manifest_sheet Name of the sheet in the excel file containing the relevant sample information. Default: the first sheet.
#' @param sample_column Name of column containing sample names that matches the values in `SampleID` in the Olink output. Default = `"SampleID"`
#' @param project_column Name of the column containing the name of the project the sample is associated with. Default = `"project"`
#' @param ... not currently used
#'
#' @returns A named list of [`"tibbles"`][`tibble::tibble`]
#'
#' @export
#'
ingest_manifest <- function(
  input,
  output_directory,
  sample_column = "SampleID",
  project_column = "Project",
  manifest_sheet = "ManifestBuilder",
  ...
) {
  .reader_func(
    input = input,
    file_extension = "excel",
    sample_column = sample_column,
    project_column = project_column,
    manifest_sheet = manifest_sheet,
    ...
  )
}


#' @title multifile_write
#' @description Write multiple files into a list and name them based on the file name without file extension
#' @details `multifile_write` encapsulates `ingest_olink_data` and `ingest_manifest` along with creating the directory
#' structure matching that as defined in the OlinkHT Standard Data Package, e.g.
#' ```
#'  .
#'  └── SDP/
#'     ├── Level_1/
#'     │   ├── plate1.parquet
#'     │   ├── plate2.parquet
#'     │   └── manifest.xlsx
#'     ├── Level_2
#'     └── code
#' ```
#'
#' @param data A named list of [`tibble::tibble`]
#' @param file_extension Output type: `"parquet"` or `"csv"` (default `"parquet"`)
#' @param proj_dir path to where the files should be saved in the standard data package output format
#'
#' @importFrom purrr imap
#' @importFrom arrow write_parquet
#' @importFrom readr write_csv
#' @importFrom stringr str_glue
#'
#' @returns None
#'
#' @export
#' @examples
multifile_write <- function(
  .data,
  file_extension = c("parquet", "csv"),
  proj_dir = NULL
) {
  file_extension <- match.arg(file_extension)

  proj_dir <- proj_dir %||% getwd()
  out_dir <- stringr::str_glue("{proj_dir}/SDP/Level_1")

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  if (file_extension == "parquet") {
    purrr::iwalk(
      .x = .data,
      .f = \(x, name) {
        arrow::write_parquet(
          x,
          sink = stringr::str_glue("{out_dir}/{name}.parquet")
        )
      }
    )
  } else if (file_extension == "csv") {
    purrr::iwalk(
      .x = .data,
      .f = \(x, name) {
        readr::write_csv(
          x = x,
          file = stringr::str_glue("{out_dir}/{name}.csv")
        )
      }
    )
  }
}


#' @title Olink_Reader
#' @description Setup project directory and import Olink data and manifest
#' @detail Setup `"input"` according to the project format described in the Olink data standards
#' document, and then read both parquet files and manfests
#' and reads them into a list format
#'
#' @param input project directory to which `"code"`, `"level 1"`, and `"level 2"` files should be
#' written and where the output from NPXExplorer and the sample manifests can be found
#' @param output path to where the standard data package directory structure should be setup and files written to
#' @inheritParams ingest_manifest manifest_sheet sample_column project_column
#'
#' @importFrom purrr map
#' @returns list with
#' * `"data"`: a list of [`tibbles`][`tibble::tibble`] with the Olink data
#' * `"manifest"`: a [`tibble`][`tibble::tibble`] with manifest information
#'
#' @export
#' @examples
#' Olink_Reader(
#'    input = "path/to/folder/with/manifest/and/parquet/files",
#'    output = "path/to/output",
#'    manifest_sheet = "manifest",
#'    sample_column = "Tube ID",
#'    project_column = "Project"
#' )
Olink_Reader <- function(
  input,
  output,
  manifest_sheet = 1,
  sample_column = NULL,
  project_column = NULL
) {
  # Defining the directory all Olink data is stored in

  # make sub-directory for SDP hiearchical file structure
  purrr::map(
    .x = c(
      root = "SDP", # root directory for SDP
      lvl1 = "SDP/Level_1", # level 1 directory for SDP
      lvl2 = "SDP/Level_2", # level 2 directory for SDP
      code = "SDP/code"
    ),
    .f = \(x) {
      if (!x %in% list.dirs(output, full.names = FALSE)) {
        dir.create(
          stringr::str_glue("{output}/{x}"),
          recursive = TRUE
        )
      } else {
        warning(stringr::str_glue(
          "An existing {x} directory was found in {output}; potentially overwriting files."
        ))
      }
    }
  )

  # Return parquet files
  list(
    data = ingest_olink_data(
      input
    ),
    manifest = ingest_manifest(
      input,
      sample_column = sample_column,
      project_column = project_column,
      manifest_sheet = manifest_sheet
    )
  )
}
