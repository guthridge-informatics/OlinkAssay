#' @title generate_performance_report
#' @description Generate a QC performance report for a Olink HT
#'  project.
#'
#' @param olink_assay
#' @param project_name title to give the report
#' @param correction_procedure description to use in
#'  the report for the batch correction procedure
#' @param output_dir path to where the report should be output
#' @param report_name name to give the output report
#' @param template_file path to Quarto qmd file to use for report.
#'
#' @returns Nothing, but writes a PDF named `report_name` to `output_dir`
#'
#' @importFrom quarto quarto_render
#' @importFrom withr local_tempfile
#'
#' @export
#' @examples
#' \dontrun{
#' generate_performance_report(
#'    olink_assay = here::here("olink_assay.RDS"),
#'    project_name = "Science Fair Project",
#'    correction_procedure = "Median correction",
#'    output_dir = "/home/user/reports",
#'    report_name = "assay_performance_report.pdf"
#'  )
#' }
generate_performance_report <- function(
  olink_assay,
  project_name = "Science Project",
  correction_procedure = "median",
  output_dir = here::here(),
  report_name = "assay_performance_report.pdf",
  template_file = NULL
) {
  template_file <- template_file %||%
    system.file(
      "qmd",
      "performance_report_template.qmd",
      package = "OlinkAssay"
    )

  tmp_path <- withr::local_tempfile()

  tmp_template_file <- glue::glue(
    "{tmp_path}{.Platform[['file.sep']]}performance_report_template.qmd"
  )
  file.copy(from = template_file, to = tmp_template_file)

  # Render report beside the template
  quarto::quarto_render(
    input = tmp_template_file,
    output_format = "pdf",
    output_file = "performance_report_template.pdf",
    execute_params = list(
      project_name = project_name,
      correction_procedure = correction_procedure,
      olink_assay_file = olink_assay
    ),
    quiet = FALSE
  )

  output_report <- glue::glue(
    "{tmp_path}{.Platform[['file.sep']]}performance_report_template.qmd"
  )

  if (file.exists(output_report)) {
    file.copy(
      from = output_report,
      to = glue::glue(
        "{sdp_directory}{.Platform[['file.sep']]}performance_report.pdf"
      )
    )
  }

  invisible(output_report)
}
