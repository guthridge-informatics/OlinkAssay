#' @title import_olink_data
#' @description Read multiple Olink output files into a list and name them based on the file name
#'
#' @param input path to where the individual Olink run data can be found
#'
#' @returns A named list of tibbles
#' @export
#' @examples
#' \dontrun{}
#' import_olink_data(input = "/path/to/data")
#'
import_olink_data <- function(input) {
  files_run <- list.files(
    path = input,
    pattern = "\\.parquet$",
    full.names = TRUE
  ) |>
    na.omit() |>
    as.character()
  if (length(files_run) < 1) {
    stop(stringr::str_glue("No {file_extension} files were found!"))
  }
  .data <- purrr::map(.x = files_run, .f = arrow::read_parquet)

  names(.data) <- stringr::str_split_i(
    string = files_run,
    pattern = .Platform$file.sep,
    i = -1
  ) |>
    stringr::str_remove(pattern = ".parquet")

  .data
}

#' @title import_manifest
#' @description Read in the manifest of samples for an Olink run
#'
#' @param manifest path to where Olink run manifest. Manifest must be in excel or comma-delimited format.
#' @param manifest_sheet Name of the sheet in the excel file containing the relevant sample information. Default: the first sheet.
#' @param sample_column Name of column containing sample names that matches the values in `SampleID` in the Olink output. Default = `"SampleID"`
#' @param project_column Name of the column containing the name of the project the sample is associated with. Default = `"project"`
#' @param additional_columns Any additional collumns from the manifest that should be carried over into later metadata
#'
#' @returns A [`"tibble"`][`tibble::tibble`] containing only the the `sample_column`, `project_column`, anything in `additional_columns`
#'
#' @export
#' @examples
#' \dontrun{}
#' import_manifest(
#'    input = "/path/to/data",
#'    sample_column = `Tube ID`,
#'    project_column = `Project`,
#'    additional_columns = c(`Sample Type`)
#' )
#'
import_manifest <- function(
  manifest,
  manifest_sheet = "ManifestBuilder",
  sample_column = "SampleID",
  project_column = "Project",
  additional_columns = NULL
) {
  if (!is.null(sample_column)) {
    sample_column <- rlang::sym(sample_column)
  }
  if (!is.null(project_column)) {
    project_column <- rlang::sym(project_column)
  }

  if (!is.null(additional_columns)) {
    loc <- tidyselect::eval_select(
      rlang::expr(additional_columns),
      data = manifest
    )
  } else {
    loc <- NULL
  }

  extension <- stringr::str_split_i(string = manifest, pattern = "\\.", i = -1)
  .data <- switch(
    EXPR = extension,
    "xlsx" = readxl::read_excel(path = manifest, sheet = manifest_sheet),
    "csv" = readr::read_csv(file = manifest)
  ) |>
    dplyr::select(
      SampleID = {{ sample_column }},
      Project = {{ project_column }},
      loc
    )

  data_name <- stringr::str_split_i(
    string = manifest,
    pattern = .Platform$file.sep,
    i = -1
  ) |>
    stringr::str_remove(pattern = glue::glue(".{extension}"))

  list(data_name = .data)
}


#' @title multifile_write
#' @description Write multiple files into a list and name them based on the file name without file extension
#' @details `multifile_write` encapsulates `import_olink_data` and `import_manifest` along with creating the directory
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


#' @title olink_reader
#' @description Setup project directory and import Olink data and manifest
#' @detail Setup `"input"` according to the project format described in the Olink data standards
#' document, and then read both parquet files and manfests
#' and reads them into a list format
#'
#' @param input project directory to which `"code"`, `"level 1"`, and `"level 2"` files should be
#' written and where the output from NPXExplorer and the sample manifests can be found
#' @param output path to where the standard data package directory structure should be setup and files written to
#' @inheritParams import_manifest manifest_sheet sample_column project_column
#'
#' @importFrom purrr map
#' @returns list with
#' * `"data"`: a list of [`tibbles`][`tibble::tibble`] with the Olink data
#' * `"manifest"`: a [`tibble`][`tibble::tibble`] with manifest information
#'
#' @export
#' @examples
#' olink_reader(
#'    input = "path/to/folder/with/manifest/and/parquet/files",
#'    output = "path/to/output",
#'    manifest_sheet = "manifest",
#'    sample_column = "Tube ID",
#'    project_column = "Project"
#' )
olink_reader <- function(
  input,
  output,
  manifest_sheet = 1,
  sample_column = NULL,
  project_column = NULL
) {
  # Defining the directory all Olink data is stored in

  setup_sdp(output)

  manifest <- list.files(input, pattern = "xlsx", full.names = TRUE)
  # Return parquet files
  list(
    data = import_olink_data(
      input
    ),
    manifest = import_manifest(
      manifest = manifest,
      sample_column = sample_column,
      project_column = project_column,
      manifest_sheet = manifest_sheet
    )
  )
}

#' @title setup_sdp
#' @description Setup a SDP hiearchical file structure
#'
#' @returns Nothing, but at `path`, produces a structure matching
#' that as defined in the OlinkHT Standard Data Package, e.g.
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
#' @export
#' @examples
#' tmp_path <- withr::local_tempfile()
#' setup_sdp(path = tmp_path)
setup_sdp <- function(path) {
  purrr::walk(
    .x = c(
      root = "SDP", # root directory for SDP
      lvl1 = "SDP/Level_1", # level 1 directory for SDP
      lvl2 = "SDP/Level_2", # level 2 directory for SDP
      code = "SDP/code"
    ),
    .f = \(x) {
      if (!x %in% list.dirs(path, full.names = FALSE)) {
        dir.create(
          stringr::str_glue("{path}/{x}"),
          recursive = TRUE
        )
      } else {
        warning(stringr::str_glue(
          "An existing {x} directory was found in {output}; potentially overwriting files."
        ))
      }
    }
  )
}
