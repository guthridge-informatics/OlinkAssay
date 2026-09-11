#' @title QCPlot
#' @description Generate plots showing Assay and Sample failure/pass rates
#' @export
setGeneric("QCPlot", function(x, ...) standardGeneric("QCPlot"))

#' @importFrom dplyr select filter row_number arrange mutate if_else
#' @importFrom ggplot2 ggplot aes geom_bar coord_polar theme_void theme scale_fill_manual geom_histogram stat_bin scale_x_continuous scale_y_continuous labs theme_classic facet_wrap vars
#' @importFrom ggrepel geom_label_repel
#' @importFrom scales percent
#'
#' @export
setMethod(
  f = "QCPlot",
  signature = "OlinkAssay",
  definition = function(x) {
    .tbl <- SampleAssayQC(x)

    plot <-
      .tbl[["AssayOlinkQC"]] |>
      dplyr::arrange(n) |>
      dplyr::mutate(
        AssayQC = factor(AssayQC, levels = c("PASS", "WARN", "FAIL")),
        prop = n / sum(n) * 100,
        ypos = cumsum(prop) - 0.5 * prop
      ) |>
      ggplot2::ggplot(ggplot2::aes(x = "", y = n, fill = AssayQC)) +
      ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
      ggplot2::coord_polar("y", start = 0) +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "none") +
      ggrepel::geom_label_repel(
        ggplot2::aes(y = ypos, label = paste0(AssayQC, "=", n)),
        size = 4,
        nudge_x = 0.1,
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(values = c("green4", "yellow2", "red")) +
      ggplot2::facet_wrap(ggplot2::vars(PlateID))

    plot_hist <-
      .tbl[["SampleOlinkQC"]] |>
      ggplot2::ggplot(mapping = ggplot2::aes(x = Frequency, fill = PlateID)) +
      ggplot2::geom_histogram(
        color = "white",
        alpha = 0.6,
        bins = 30
      ) +
      ggplot2::stat_bin(
        bins = 30,
        geom = "text",
        mapping = ggplot2::aes(
          label = dplyr::if_else(
            condition = ggplot2::after_stat(count) == 0,
            true = "",
            false = as.character(ggplot2::after_stat(count))
          )
        ),
        vjust = -0.3,
        size = 4
      ) +
      ggplot2::scale_x_continuous(
        labels = scales::percent,
        limits = c(-0.05, 1)
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.2))
      ) +
      ggplot2::labs(
        x = "% analytes failed",
        y = "Number of samples"
      ) +
      ggplot2::theme_classic() +
      ggplot2::facet_wrap(ggplot2::vars(PlateID))

    cowplot::plot_grid(plotlist = list(plot, plot_hist), nrow = 2)
  }
)


#' @title batchCorrectionUMAP
#' @description Prepare a UMAP graph from level 2 data to examine the success or failure of batch correction
#'
#' @param exclude_high_variance_assays should high variance assays be included when calculating the UMAP?
#' Default: FALSE
#' @param title title to add to plot
#'
#' @returns A ggproto object (e.g. the plot object)
#'
#' @export
#' @examples
setGeneric("batchCorrectionUMAP", function(x, ...) {
  standardGeneric("batchCorrectionUMAP")
})

#' @rdname batchCorrectionUMAP
#' @param x OlinkAssay object containing batch-corrected data
#' @param assay_name assay to use for expression data. One of `"ExtNPX"`, `"ExtNPX_Corrected"`, `"LogProtExp_Raw"`
#'
#' @importFrom glue glue
#' @importFrom SummarizedExperiment assayNames asssay
#' @importFrom dplyr select mutate across if_else rename left_join join_by
#' @importFrom tidyselect where
#' @importFrom tibble as_tibble
#' @importFrom ggplot2 ggplot aes geom_point theme theme_minimal element_text element_blank ggtitle
#' @importFrom uwot umap
#' @export
setMethod(
  f = "batchCorrectionUMAP",
  signature = "OlinkAssay",
  definition = function(
    x,
    assay_name = c("ExtNPX", "ExtNPX_Corrected", "LogProtExp_Raw"),
    exclude_high_variance_assays = FALSE,
    title = NULL
  ) {
    assay_name <- match.arg(assay_name)
    if (!assay_name %in% SummarizedExperiment::assayNames(x)) {
      stop(glue::glue(
        "{.assay} was not found in the object. Available assays include {glue::glue_collapse(assayNames(x), sep=', ', last = ', and ')}"
      ))
    }
    # .col <- rlang::sym(.col)

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
      SummarizedExperiment::assay(x, assay_name)[
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
  }
)
