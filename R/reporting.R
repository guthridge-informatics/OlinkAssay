# TODO: replace rmarkdown with quarto

#' @title generate_performance_report
#'
#' @param project_name
#' @param correction_procedure
#' @param sdp_directory
#' @param report_name
#' @param template_dir
#'
#' @returns
#'
#' @importFrom rmarkdown render
#'
#' @export
#' @examples
generate_performance_report <- function(
  project_name,
  correction_procedure,
  sdp_directory,
  output_dir,
  report_name = "assay_performance_report.pdf"
) {
  # Construct commonly used paths

  template_file <- system.file(
    "rmd",
    "performance_report_template.Rmd",
    package = "olinkqc"
  )

  output_pdf <- file.path(
    sdp_directory,
    report_name
  )

  level_1_path <- file.path(
    sdp_directory,
    "Level_1"
  )

  level_2_path <- file.path(
    sdp_directory,
    "Level_2"
  )

  # Render report beside the template
  rmarkdown::render(
    input = template_file,
    output_format = "pdf_document",
    output_dir = output_dir,
    params = list(
      project_name = project_name,
      correction_procedure = correction_procedure,
      level_1_file_path = level_1_path,
      level_2_file_path = level_2_path
    ),
    quiet = TRUE
  )

  # Accessing the rendered pdf in the temporary storage location
  rendered_pdf <- file.path(
    output_dir,
    "performance_report_template.pdf"
  )

  # Copy finished report to SDP directory
  # fs::file_copy(
  #   path = rendered_pdf,
  #   new_path = output_pdf,
  #   overwrite = TRUE
  # )

  invisible(output_pdf)
}
