# optimize the UMAP visualization for any omics

#' @title vst_to_pca
#' @description Use variance stabilized data to detect variables/analytes with above expected variances
#'
#' @param .data A sample-by-assay [`base::data.frame`] or [`tibble::tibble`] with variance-stablized data
#' @param exclude Assays to ignore when calculating the PCA
#'
#' @returns list with:
#' 1. `data_vst` = [`tibble::tibble`] with assay mean, variance, predicted variance, and difference score above the predicted means
#' 2. `plot` = [`ggplot2::ggplot`] of mean vs variance with the loess function overlayed
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth theme_classic
#' @importFrom dplyr select across summarise
#' @importFrom tidyselect where
#' @importFrom stats loess
#'
#' @export
#' @examples
vst_to_pca <- function(.data, exclude = NULL) {
  # data excludes the categorical columns that will not be used for the vst process
  .data <- dplyr::select(.data, -tidyselect(all_of(exclude)))

  # calculate assay mean and variance
  if (tibble::is_tibble(.data)) {
    # mean <- apply(X = data, MARGIN = 2, FUN = \(x) mean(unlist(x)))
    data_vst <-
      .data |>
      tidyr::pivot_longer(
        cols = tidyselect::where(is.numeric),
        names_to = "Assay"
      ) |>
      dplyr::group_by(Assay) |>
      dplyr::summarise(Mean = mean(value), Variance = var(value))
  } else {
    data_vst <- tibble::tibble(
      Assay = colnames(.data),
      Mean = colMeans(x = .data),
      Variance = colMeans(x = .data)
    )
  }

  # plot the scatter plots of the mean and variance of each protein
  fit_lowess <- stats::loess(Variance ~ Mean, data_vst, span = 0.2)

  # graph the loess function in the plot of mean vs variance
  plot <-
    ggplot2::ggplot(
      data = data_vst,
      mapping = ggplot2::aes(x = Mean, y = Variance)
    ) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(
      method = loess,
      formula = y ~ x,
      lty = 2,
      method.args = list(span = 0.2)
    ) +
    ggplot2::theme_classic()

  # add predicted variance
  # data_vst[["PredictedVar"]] <- predict(fit_lowess, data_vst$Mean)

  # add difference score above the predicted means
  # data_vst[["Diff"]] <- data_vst[["Variance"]] - data_vst[["PredictedVar"]]

  data_vst <-
    data_vst |>
    dplyr::mutate(
      PredictedVar = predict(fit_lowess, Mean),
      Diff = Variance - PredictedVar
    )

  list(
    data_vst = data_vst,
    plot = plot
  )
}

#' @title pca_to_umap
#' @description Automatically find the best PCA dimensional cut-off for the subsequent UMAP
#'
#' @param .data
#'
#' @returns
#'
#' @importFrom stats prcomp
#' @importFrom dplyr mutate
#' @importFrom pathviewr find_curve_elbow
#'
#' @export
#' @examples
pca_to_umap <- function(.data) {
  # calculate the principal components
  pca <- prcomp(stats::na.omit(.data), scale. = TRUE, center = TRUE)
  # get PCA importance from the principal component analysis
  res_pca <- summary(pca)$importance |>
    t() |>
    data.frame() |>
    dplyr::mutate(PCs = c(seq_len(nrow(.))))
  n_pca <- pathviewr::find_curve_elbow(
    data_frame = res_pca[, c("PCs", "Proportion.of.Variance")],
    plot_curve = TRUE
  )

  pca[["x"]][, c(1:n_pca)] |> data.frame()
}

# There has to be a way of using one function below and then just swapping the argument that is matched, but damned if I can figure out {rlang}
#' @title optimize_n_neighbor
#' @description Optimization function for n_neighbors by setting the spread = 10, min_dist = 0.1 as default by can be changed as needed
#'
#' @param .data
#' @param groups
#' @param spread
#' @param min_dist
#' @param min
#' @param max
#' @param step
#'
#' @returns
#'
#' @importFrom purrr map
#' @importFrom uwot umap
#' @importFrom tibble as_tibble
#' @importFrom dplyr rename
#' @importFrom ggplot2 ggplot aes geom_point theme_classic ggtitle theme element_text
#' @importFrom stringr str_glue
#' @importFrom cowplot plot_grid
#'
#' @export
#' @examples
optimize_n_neighbor <- function(
  .data,
  groups = NULL,
  spread = 10,
  min_dist = 0.1,
  min = 5,
  max = 16,
  step = 1
) {
  if (is.null(groups)) {
    groups <- "black"
  }
  plots <- purrr::map(
    .x = seq(min, max, step),
    .f = \(n_neighbors) {
      umap_cyto <- uwot::umap(
        .data,
        spread = spread,
        min_dist = min_dist,
        n_neighbors = n_neighbors,
        seed = 123
      )
      dat_umap <- umap_cyto |>
        tibble::as_tibble() |>
        dplyr::rename("UMAP1" = "V1", "UMAP2" = "V2")

      plot_0 <- ggplot2::ggplot(
        data = dat_umap,
        mapping = ggplot2::aes(
          x = UMAP1,
          y = UMAP2
        )
      ) +
        ggplot2::geom_point(
          ggplot2::aes(color = groups),
          size = 4
        ) +
        ggplot2::theme_classic() +
        ggplot2::ggtitle(stringr::str_glue("n_neighbor = {n_neighbors}")) +
        ggplot2::theme(
          legend.text = ggplot2::element_text(size = 12),
          axis.text = ggplot2::element_text(size = 12)
        )
    }
  )
  cowplot::plot_grid(plotlist = plots, ncol = 4)
}

#' @title optimize_spread
#' @description Optimization function for spread after n_neighbor parameter has been fixed, by default this function uses min_dist of 0.1
#'
#' @param .data
#' @param groups
#' @param n_neighbor
#' @param min_dist
#' @param min
#' @param max
#' @param step
#'
#' @returns
#'
#' @importFrom purrr map
#' @importFrom uwot umap
#' @importFrom tibble as_tibble
#' @importFrom dplyr rename
#' @importFrom ggplot2 ggplot aes geom_point theme_classic ggtitle theme element_text
#' @importFrom stringr str_glue
#' @importFrom cowplot plot_grid
#'
#' @export
#' @examples
optimize_spread <- function(
  .data,
  groups,
  n_neighbor,
  min_dist = 0.1,
  min = 1,
  max = 15,
  step = 1
) {
  if (is.null(groups)) {
    groups <- "black"
  }
  plots <- purrr::map(
    .x = seq(min, max, step),
    .f = \(spread) {
      umap_cyto <- uwot::umap(
        .data,
        spread = spread,
        min_dist = min_dist,
        n_neighbors = n_neighbors,
        seed = 123
      )
      dat_umap <- umap_cyto |>
        tibble::as_tibble() |>
        dplyr::rename("UMAP1" = "V1", "UMAP2" = "V2")

      plot_0 <- ggplot2::ggplot(
        data = dat_umap,
        mapping = ggplot2::aes(
          x = UMAP1,
          y = UMAP2
        )
      ) +
        ggplot2::geom_point(
          ggplot2::aes(color = groups),
          size = 4
        ) +
        ggplot2::theme_classic() +
        ggplot2::ggtitle(stringr::str_glue("n_neighbor = {n_neighbors}")) +
        ggplot2::theme(
          legend.text = ggplot2::element_text(size = 12),
          axis.text = ggplot2::element_text(size = 12)
        )
    }
  )
  cowplot::plot_grid(plotlist = plots, ncol = 4)
}

#' @title optimize_min_dist
#' @description The final optimization function that is to fix the min_dist after n_neighbor and spread parameters have been fixed
#'
#' @param .data
#' @param groups
#' @param spread
#' @param n_neighbor
#' @param min
#' @param max
#' @param step
#'
#' @returns
#'
#' @importFrom purrr map
#' @importFrom uwot umap
#' @importFrom tibble as_tibble
#' @importFrom dplyr rename
#' @importFrom ggplot2 ggplot aes geom_point theme_classic ggtitle theme element_text
#' @importFrom stringr str_glue
#' @importFrom cowplot plot_grid
#'
#' @export
#' @examples
optimize_min_dist <- function(
  .data,
  groups,
  spread,
  n_neighbor,
  min = 0.01,
  max = 0.5,
  step = 0.03
) {
  if (is.null(groups)) {
    groups <- "black"
  }
  plots <- purrr::map(
    .x = seq(min, max, step),
    .f = \(min_dist) {
      umap_cyto <- uwot::umap(
        .data,
        spread = spread,
        min_dist = min_dist,
        n_neighbors = n_neighbors,
        seed = 123
      )
      dat_umap <- umap_cyto |>
        tibble::as_tibble() |>
        dplyr::rename("UMAP1" = "V1", "UMAP2" = "V2")

      plot_0 <- ggplot2::ggplot(
        .data = dat_umap,
        mapping = ggplot2::aes(
          x = UMAP1,
          y = UMAP2
        )
      ) +
        ggplot2::geom_point(
          ggplot2::aes(color = groups),
          size = 4
        ) +
        ggplot2::theme_classic() +
        ggplot2::ggtitle(stringr::str_glue("n_neighbor = {n_neighbors}")) +
        ggplot2::theme(
          legend.text = ggplot2::element_text(size = 12),
          axis.text = ggplot2::element_text(size = 12)
        )
    }
  )
  cowplot::plot_grid(plotlist = plots, ncol = 4)
}
