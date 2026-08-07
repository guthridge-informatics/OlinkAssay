#' @title normalization_check
#'
#' @param data_corrected
#' @param pt.size
#'
#' @returns
#'
#' @importFrom dplyr filter rename mutate select
#' @importFrom tidyr pivot_wider
#' @importFrom ggplot2 ggtitle
#' @export
#' @examples
normalization_check <- function(data_corrected, pt.size = 0.5) {
  data_corrected_combined <-
    data_corrected |>
    dplyr::filter(AssayType == "assay" & SampleType == "SAMPLE") |>
    dplyr::rename("RawExtNPX" = "ExtNPX") |>
    dplyr::mutate(ProteinID = paste0(Assay, "_", OlinkID)) |>
    dplyr::select(
      SampleID,
      PlateID,
      ProteinID,
      RawExtNPX,
      ExtNPX_Corrected,
      LogProtExp_Raw
    ) |>
    na.omit() |>
    tidyr::pivot_wider(
      names_from = ProteinID,
      values_from = c(RawExtNPX, ExtNPX_Corrected, LogProtExp_Raw)
    )

  rawextnpx_data <- data_corrected_combined |>
    dplyr::select(contains("RawExtNPX"))

  p1 <- rawextnpx_data |>
    na.omit() |>
    UMAP_groups(
      groups = na.omit(data_corrected_combined)[["PlateID"]],
      n_neighbors = 30,
      pt.size = pt.size
    ) +
    ggplot2::ggtitle(paste0(
      "Raw ExtNXP - ",
      length(colnames(rawextnpx_data)),
      " Proteins Visualized"
    ))

  corrextnpx_data <- data_corrected_combined |>
    dplyr::select(contains("ExtNPX_Corrected"))

  p2 <- corrextnpx_data |>
    na.omit() |>
    UMAP_groups(
      groups = na.omit(data_corrected_combined)[["PlateID"]],
      n_neighbors = 30,
      pt.size = pt.size
    ) +
    ggplot2::ggtitle(
      paste0(
        "Batch-corrected ExtNPX - ",
        length(colnames(corrextnpx_data)),
        " Proteins Visualized"
      )
    )

  logprotexp_data <- data_corrected_combined |>
    dplyr::select(dplyr::contains("LogProtExp_Raw"))

  p3 <- logprotexp_data |>
    na.omit() |>
    UMAP_groups(
      groups = na.omit(data_corrected_combined)[["PlateID"]],
      pt.size = pt.size,
      n_neighbors = 30
    ) +
    ggplot2::ggtitle(
      paste0(
        "Batch-corrected LogProtExp - ",
        length(colnames(logprotexp_data)),
        " Proteins Visualized"
      )
    )

  list(
    plot1 = p1,
    plot2 = p2,
    plot3 = p3
  )
}

#' @title
#' @description graph the clean reworked umap
#'
#' @param data
#' @param k
#' @param eps
#' @param minPts
#' @param arrow_size
#' @param pt.size
#' @param arrowtip_size
#' @param cols
#' @param label
#' @param label.size
#' @param spread
#' @param min_dist
#' @param n_neighbors
#'
#' @returns
#'
#' @importFrom uwot umap
#' @importFrom tibble as_tibble
#' @importFrom dplyr rename mutate
#' @importFrom dbscan sNNclust
#' @importFrom ggplot2 aes element_blank element_text geom_point geom_segment geom_text get_layer_data ggplot ggtitle guides guide_legend theme theme_void
#' @importFrom stringr str_glue
#' @importFrom grid arrow unit
#'
#' @export
#' @examples
UMAP_snn <- function(
  data,
  k = NA,
  eps = 7,
  minPts = 10,
  arrow_size = 0.1,
  pt.size = 4,
  arrowtip_size = 2,
  cols = NULL,
  label = FALSE,
  label.size = 15,
  spread = 5,
  min_dist = 0.25,
  n_neighbors = NA
) {
  if (is.na(n_neighbors)) {
    n_neighbors <- nrow(data) / 10
  }
  dat_umap <- uwot::umap(
    X = data,
    spread = spread,
    min_dist = min_dist,
    n_neighbors = n_neighbors,
    seed = 123
  ) |>
    tibble::as_tibble() |>
    dplyr::rename("UMAP1" = "V1", "UMAP2" = "V2")

  # sNN clustering the UMAP coordinates
  if (is.na(k)) {
    k <- nrow(dat_umap) / 10
  }
  groups <- as.factor(
    dbscan::sNNclust(dat_umap, k = k, eps = eps, minPts = minPts)[["cluster"]]
  )

  # graph the umap plot
  p <- ggplot2::ggplot(
    data = dat_umap,
    mapping = ggplot2::aes(x = UMAP1, y = UMAP2)
  ) +
    ggplot2::geom_point(ggplot2::aes(color = groups), size = 4) +
    ggplot2::theme_void()

  layer_data <- ggplot2::get_layer_data(p)
  y_min <- min(layer_data[["y"]])
  y_max <- max(layer_data[["y"]])
  x_min <- min(layer_data[["x"]])
  x_max <- max(layer_data[["x"]])

  p <- p +
    ggplot2::guides(
      fill = ggplot2::guide_legend(title = stringr::str_glue("k = {k}"))
    ) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::geom_segment(
      x = x_min,
      y = y_min,
      xend = x_min,
      yend = y_min + (y_max - y_min) * arrow_size,
      size = 0.8,
      arrow = grid::arrow(
        length = grid::unit(arrowtip_size, "mm"),
        type = "closed"
      )
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = x_min + 0.2, y = y_min, label = "UMAP1"),
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    ggplot2::geom_segment(
      x = x_min,
      y = y_min,
      xend = x_min + (x_max - x_min) * arrow_size,
      yend = y_min,
      size = 0.8,
      arrow = grid::arrow(
        length = grid::unit(arrowtip_size, "mm"),
        type = "closed"
      )
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = x_min - 0.6, y = y_min, label = "UMAP2"),
      angle = 90,
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        ncol = 1,
        label.theme = ggplot2::element_text(face = "bold", size = label.size),
        override.aes = list(size = 6)
      )
    ) +
    ggplot2::ggtitle(label = "")

  list(
    plot = p,
    data_clust = dplyr::mutate(dat_umap, Cluster = groups)
  )
}

#' @title UMAP_groups
#' @description graph the clean reworked umap
#'
#' @param data
#' @param groups
#' @param eps
#' @param minPts
#' @param arrow_size
#' @param pt.size
#' @param arrowtip_size
#' @param cols
#' @param label
#' @param label.size
#' @param spread
#' @param min_dist
#' @param n_neighbors
#'
#' @returns
#'
#' @importFrom uwot umap
#' @importFrom tibble as_tibble
#' @importFrom dplyr rename
#' @importFrom ggplot2 aes element_blank element_text geom_point geom_segment geom_text get_layer_data ggplot ggtitle guides guide_legend theme theme_void
#' @importFrom stringr str_glue
#' @importFrom grid arrow unit
#'
#' @export
#' @examples
UMAP_groups <- function(
  data,
  groups,
  eps = 7,
  minPts = 10,
  arrow_size = 0.1,
  pt.size = 0.5,
  arrowtip_size = 2,
  cols = NULL,
  label = FALSE,
  label.size = 15,
  spread = 5,
  min_dist = 0.25,
  n_neighbors = NA
) {
  if (is.na(n_neighbors)) {
    n_neighbors <- nrow(data) / 10
  }
  dat_umap <- uwot::umap(
    data,
    spread = spread,
    min_dist = min_dist,
    n_neighbors = n_neighbors,
    seed = 123
  ) |>
    tibble::as_tibble() |>
    dplyr::rename("UMAP1" = "V1", "UMAP2" = "V2")

  # graph the umap plot
  p <- ggplot2::ggplot(
    data = dat_umap,
    mapping = ggplot2::aes(x = UMAP1, y = UMAP2)
  ) +
    ggplot2::geom_point(ggplot2::aes(color = groups), size = 4) +
    ggplot2::theme_void()
  layer_data <- ggplot2::get_layer_data(p)
  y_min <- min(layer_data[["y"]])
  y_max <- max(layer_data[["y"]])
  x_min <- min(layer_data[["x"]])
  x_max <- max(layer_data[["x"]])

  p +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::geom_segment(
      x = x_min,
      y = y_min,
      xend = x_min,
      yend = y_min + (y_max - y_min) * arrow_size,
      size = 0.8,
      arrow = grid::arrow(
        length = grid::unit(arrowtip_size, "mm"),
        type = "closed"
      )
    ) +
    ggplot2::geom_text(
      mapping = ggplot2::aes(
        x = x_min + (x_max - x_min) * 0.01,
        y = y_min,
        label = "UMAP1"
      ),
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    ggplot2::geom_segment(
      x = x_min,
      y = y_min,
      xend = x_min + (x_max - x_min) * arrow_size,
      yend = y_min,
      size = 0.8,
      arrow = grid::arrow(
        length = grid::unit(arrowtip_size, "mm"),
        type = "closed"
      )
    ) +
    ggplot2::geom_text(
      mapping = ggplot2::aes(
        x = x_min - (x_max - x_min) * 0.02,
        y = y_min + (x_max - x_min) * 0.01,
        label = "UMAP2"
      ),
      angle = 90,
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        ncol = 1,
        label.theme = ggplot2::element_text(face = "bold", size = label.size),
        override.aes = list(size = 6)
      )
    ) +
    ggplot2::ggtitle(label = "")
}

#' @title graph_UMAP
#' @description graphing pretty umap graph
#'
#' @param data_umap
#' @param groups
#' @param arrow_size
#' @param pt.size
#' @param arrowtip_size
#' @param cols
#' @param label
#' @param label.size
#'
#' @returns
#'
#' @importFrom ggplot2 aes element_blank element_text geom_point geom_segment geom_text get_layer_data ggplot ggtitle guides guide_legend theme theme_void
#' @importFrom grid arrow unit
#'
#' @export
#' @examples
graph_UMAP <- function(
  data_umap,
  groups,
  arrow_size = 0.1,
  pt.size = 0.5,
  arrowtip_size = 2,
  cols = NULL,
  label = FALSE,
  label.size = 15
) {
  # graph the umap plot
  p <- ggplot2::ggplot(
    data = data_umap,
    mapping = ggplot2::aes(x = UMAP1, y = UMAP2)
  ) +
    ggplot2::geom_point(mapping = aes(color = as.factor(groups)), size = 0.5) +
    ggplot2::theme_void()
  # y-range
  yrange <- layer_scales(p)$y$range$range
  # x-range
  xrange <- layer_scales(p)$x$range$range
  p +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::geom_segment(
      x = x_min,
      y = y_min,
      xend = x_min,
      yend = y_min + (y_max - y_min) * arrow_size,
      size = 0.8,
      arrow = arrow(length = unit(arrowtip_size, "mm"), type = "closed")
    ) +
    ggplot2::geom_text(
      mapping = ggplot2::aes(
        x = x_min + (x_max - x_min) * 0.01,
        y = y_min,
        label = "UMAP1"
      ),
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    ggplot2::geom_segment(
      x = x_min,
      y = y_min,
      xend = x_min + (x_max - x_min) * arrow_size,
      yend = y_min,
      size = 0.8,
      arrow = grid::arrow(
        length = grid::unit(arrowtip_size, "mm"),
        type = "closed"
      )
    ) +
    ggplot2::geom_text(
      mapping = ggplot2::aes(
        x = x_min - (x_max - x_min) * 0.02,
        y = y_min + (x_max - x_min) * 0.01,
        label = "UMAP2"
      ),
      angle = 90,
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        ncol = 1,
        label.theme = ggplot2::element_text(face = "bold", size = label.size),
        override.aes = list(size = 6)
      )
    ) +
    ggplot2::ggtitle(label = "")
}

#' @title plot_Qc
#' @description graph the QC performance of each plate
#'
#' @param data_i
#'
#' @returns
#'
#' @importFrom dplyr select filter row_number arrange mutate if_else
#' @importFrom ggplot2 ggplot aes geom_bar coord_polar theme_void theme scale_fill_manual geom_histogram stat_bin scale_x_continuous scale_y_continuous labs theme_classic facet_wrap vars
#' @importFrom ggrepel geom_label_repel
#'
#' @export
#' @examples
plot_qc <- function(data_i) {
  plateID1 <- data_i[["AssayOlinkQC"]] |>
    dplyr::select(PlateID) |>
    dplyr::filter(dplyr::row_number() == 1) |>
    unique()

  plot <-
    data_i[["AssayOlinkQC"]] |>
    dplyr::filter(PlateID == plateID1[["PlateID"]]) |>
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
    ggplot2::geom_label_repel(
      ggplot2::aes(y = ypos, label = paste0(AssayQC, "=", n)),
      size = 4,
      nudge_x = 0.1,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = c("green4", "yellow2", "red"))

  plot_hist <-
    data_i[["SampleOlinkQC"]] |>
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
          false = ggplot2::after_stat(count)
        )
      ),
      vjust = -0.3,
      size = 4
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::percent,
      limits = c(-0.05, 1)
    ) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    ggplot2::labs(
      x = "% analytes failed",
      y = "Number of samples"
    ) +
    ggplot2::theme_classic() +
    ggplot2::facet_wrap(ggplot2::vars(PlateID))

  cowplot::plot_grid(plotlist = list(plot, plot_hist), nrow = 2)
}
