entry_point <- function(
    manifest,
    manifest_sheet,
    project_directory,
    output_directory = NULL,
    projects_to_include,
    project_name,
    correction_procedure,
    report_name = NULL,
    kitlot,
    pclot,
    sclot,
    nclot,
    batch_correction_method = c("median", "global")
) {
    batch_correction_method <- match.arg(batch_correction_method)

    output_directory <- output_directory %||% project_directory

    if (!"bridge" %in% projects_to_include) {
        projects_to_include <- c("bridge", projects_to_include)
    }

    manifest_information <- readxl::read_excel(
        path = manifest,
        sheet = manifest_sheet
    )

    # Gathering Files and creating the SDP
    combined_files <- olink_reader(
        input = project_directory,
        output = output_directory,
        manifest_sheet = manifest_sheet
    ) # one *really* should pass arguments for the `sample_column`, `project_column`, and `manifest_sheet` here

    level1_data <- olink_lvl1(
        olink_files = combined_files,
        # Project names derived from manifests
        proj_names = projects_to_include, # I guess we could include "Bridge" even if it isn't an argument?
        # Lot Information
        KitLot_info = kitlot,
        PCLot_info = pclot,
        SCLot_info = sclot,
        NCLot_info = nclot
    )

    # Performing batch correction using the median Plate Control Correction Procedure
    batch_corrected_data <- batch_correction(
        level1_data[["data"]],
        method = batch_correction_method
    )

    # Preparing the data for Level 2 QC
    level_2_prep <- purrr::map(.x = batch_corrected_data, .f = olink_lvl2_prep)

    # Performing Level 2 QC
    level_2_data <- olink_lvl2(level_2_prep)

    if (is.na(correction_procedure)) {
        report_correction_method <- switch(
            EXPR = batch_corrected_method,
            median = "Plate Control Batch Correction",
            global = "Global Median Batch Correction"
        )
    }

    report_name <- report_name %||%
        glue::glue(
            "{glue::glue_collapse(projects[which(projects != 'Bridge')],sep=', ')} Olink HT Performance Report"
        )
    # Generation, Assay Performance Report
    generate_performance_report(
        project_name = project_name,
        correction_procedure = report_correction_method,
        sdp_directory = glue::glue(
            "{project_directory}{.Platform[['file.sep']]}SCP"
        ),
        report_name = report_name
    )

    n_analytes <- level_2_data$Assay |>
        unique() |>
        length()

    # Preprocessing for manifest generation
    level1_data[["data"]] |>
        # Converting level 1 filename to the Assay_Filename column
        dplyr::bind_rows(.id = "Assay_Filename") |>
        # Filtering to only Technical Bridges and Samples
        dplyr::filter(SampleType == "SAMPLE") |>
        # Counting SampleQC entries per columnn
        dplyr::count(
            Assay_Filename,
            SampleID,
            Project,
            SampleQC
        ) |>
        # Entering SampleQC field Counts as separate columns
        tidyr::pivot_wider(
            names_from = SampleQC,
            values_from = n,
            values_fill = 0
        ) |>
        # Counting the percentage of Fail and Warn analytes per Sample
        mutate(
            FAIL = dplyr::coalesce(FAIL, 0),
            WARN = dplyr::coalesce(WARN, 0),
            Analyte_Failure_Pct = 100 * FAIL / n_analytes,
            Analyte_Warn_Pct = 100 * WARN / n_analytes
        ) |>
        # Joining the manifest information from the original received manifest
        dplyr::left_join(
            manifest_information,
            by = dplyr::join_by(SampleID == `Tube ID`)
        ) |>
        # Fields not included in the original manifest
        dplyr::mutate(
            Same_Day_Processing = NA,
            Serum_Hemolysis = NA,
            Serum_Lipemia = NA,
            Processing_Comments = NA
        ) |>
        # Selecting and Ordering the Columns Prior to Writing the Manifest
        dplyr::select(
            SubjectRef,
            VisitRef,
            `Sample Alias`,
            `Sample Type`,
            Same_Day_Processing,
            Serum_Hemolysis,
            Serum_Lipemia,
            Processing_Comments,
            Assay_Filename,
            SampleID,
            Analyte_Failure_Pct,
            Analyte_Warn_Pct,
            Project
        ) |>
        # Reformatting Manifest Column Names
        dplyr::rename(
            "Assay_Sample_ID" = "SampleID",
            "Sample_Alias" = "Sample Alias",
            "Sample_Type" = "Sample Type"
        ) |>
        # Writing a preliminary excel template to be cleaned
        # according to the format of the manifest_template.xlsx file in the resources folder
        readr::write_csv(glue::glue(
            "{output_directory}{.Platform[['file.sep']]}manifest_raw.csv"
        ))
}
