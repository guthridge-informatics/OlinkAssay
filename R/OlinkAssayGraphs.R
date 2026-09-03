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
