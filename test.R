x <- readFromDisk(
  npx_file = here::here(
    "tests",
    "initial_data",
    "FS19030710_Extended_NPX_2026-06-22.parquet"
  ),
  metadata_file = here::here(
    "tests",
    "initial_data",
    "2026-05-15 James ABC-PreEVADE_Peds Serum_Olink (C.Guthridge)_Manifest.xlsx"
  ),
  metadata_sheet = "ManifestBuilder",
  sample_column = "Tube ID",
  verbose = TRUE
)
y <- readFromDisk(
  npx_file = here::here(
    "tests",
    "initial_data",
    "FS19030751_Extended_NPX_2026-06-22.parquet"
  ),
  metadata_file = here::here(
    "tests",
    "initial_data",
    "2026-05-15 James ABC-PreEVADE_Peds Serum_Olink (C.Guthridge)_Manifest.xlsx"
  ),
  metadata_sheet = "ManifestBuilder",
  sample_column = "Tube ID",
  verbose = TRUE
)


x <- oa
.assay <- "ExtNPX"
exclude_high_variance_assays <- FALSE
title <- "It works!"

if (exclude_high_variance_assays) {
  included_assays <- rownames(rowData(x))[which(
    rowData(x)[["AssayType"]] == "assay" &
      rowData(x)[["high_var_assay"]] == "Pass"
  )]
} else {
  included_assays <- rownames(rowData(x))[which(
    rowData(x)[["AssayType"]] == "assay"
  )]
}

included_samples <-
  rownames(colData(x)[which(colData(x)[["Project"]] != "Bridge"), ])

umap_intermediate <-
  assay(x, .assay)[
    included_assays,
    included_samples
  ]

sample_md <-
  colData(oa) |>
  tibble::as_tibble(rownames = "SampleID") |>
  dplyr::select(
    SampleID,
    PlateID,
    Project
  )

p <- umap_intermediate |>
  t() |>
  as.data.frame() |>
  dplyr::mutate(
    dplyr::across(
      .cols = tidyselect::where(is.numeric),
      .fns = \(x) {
        dplyr::if_else(
          is.nan(x),
          true = 0,
          false = x
        )
      }
    )
  ) |>
  uwot::umap(
    scale = TRUE,
    min_dist = 0.4,
    seed = 825
  ) |>
  tibble::as_tibble(
    rownames = "SampleID",
    .name_repair = "universal_quiet"
  ) |>
  dplyr::rename(
    UMAP1 = `...1`,
    UMAP2 = `...2`
  ) |>
  dplyr::left_join(
    y = sample_md,
    by = dplyr::join_by(SampleID == SampleID)
  ) |>
  ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = UMAP1,
      y = UMAP2,
      color = PlateID
    )
  ) +
  ggplot2::geom_point(size = 3) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.text = ggplot2::element_text(face = "bold", size = 12),
    legend.title = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = 12)
  )

if (!is.null(title)) {
  p + ggplot2::ggtitle(title)
} else {
  p
}

quarto::quarto_render(
  input = template_file,
  output_format = "pdf",
  output_file = "performance_report_template.pdf",
  execute_params = list(
    project_name = project_name,
    correction_procedure = correction_procedure,
    olink_assay_file = olink_assay
  ),
  quiet = FALSE
)
