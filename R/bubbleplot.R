#' bubbleplot
#'
#' Display percentage of cells expressing and levels for a set of genes.
#'
#' @param seuratObj Seurat object
#' @param genes_plot Either a list of genes or a data frame of annotated genes
#'   to display (see 'annotated_gene_list'). Note: Gene and protein names
#'   may be converted to the proper gene name automagically by
#'   HGNChelper::checkGeneSymbols if 'translate_gene_names' is TRUE.
#'   Genes not appearing in the dataset are skipped.
#' @param filter_exp_pct Display only genes that are expressed above this
#'   fraction of cells in at least one group. Default: NULL
#' @param filter_exp_pct_thresh Threshold for expression fraction. Default: 0
#' @param filter_exp_level Display only genes that are expressed above this
#'   level in at least one group. Default: 0
#' @param group_by Variable by which to group cells.  Can be cluster identities
#'   or any column from the meta.data slot. Default: ident
#' @param x_lab_size Font size for the x-axis labels. Default: 9
#' @param y_lab_size Font size for the y-axis labels. Default: 9
#' @param x_lab_rot_angle Angle to rotate the x-axis labels. Default: 45°
#' @param order_genes Should the genes by displayed in alphabetical order? Default: FALSE.
#' @param cluster_x Arrange the x-axis variables using hierarchical clustering.
#'   Default: TRUE
#' @param cluster_y Arrange the y-axis variables using hierarchical clustering.
#'   Default: FALSE
#' @param colors_use Color palette to use to display expression levels.
#'   Default: "Reds"
#' @param translate_gene_names Should gene names be checked and translated to
#'   their correct HUGO Gene Nomenclature Committee name? Default: FALSE
#' @param annotated_gene_list Is the gene list annotated?  If so, the genes will be
#'   faceted using their annontations.  Requires that 'genes_plot' is a two column
#'   with the annotations in a column named 'annotations' and genes in a column
#'   named 'genes'. Default: FALSE
#' @param do_return Return a ggplot2 object instead of displaying
#'
#' @import ggplot2
#' @import magrittr
#' @importFrom dplyr group_by summarize mutate ungroup
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom Seurat FetchData AverageExpression SetAllIdent
#' @importFrom tidyr gather
#' @importFrom stats hclust dist as.dendrogram order.dendrogram
#' @importFrom compositions normalize
#' @importFrom gtools mixedorder
#' @importFrom HGNChelper checkGeneSymbols
#' @importFrom glue glue
#'
#' @return if isTRUE(do_return), a ggplot2 object
#' @export
#'
#' @examples BubblePlot(seuratObj = obj, genes_plot = c("IFIT1","IFITM1","IFITM3"), group_by = "treatment")
bubbleplot <- function(seuratObj,
                       genes_plot,
                       use_scaled = FALSE,
                       filter_exp_pct = NULL,
                       filter_exp_pct_thresh = 0,
                       filter_exp_level = 0,
                       group_by = "ident",
                       x_lab_size = 9,
                       y_lab_size = 9,
                       x_axis_title = "Genes",
                       y_axis_title = "Grouping",
                       pct_legend_title = "Percent group expressing",
                       scale_legend_title = "Average scaled expression",
                       x_lab_rot_angle = 45,
                       order_genes = FALSE,
                       cluster_x = TRUE,
                       cluster_y = FALSE,
                       colors_use = NULL,
                       translate_gene_names = FALSE,
                       annotated_gene_list = FALSE,
                       do_return = FALSE) {
  if (annotated_gene_list) {
    genes_list <- genes_plot
    genes_plot <- genes_list$genes
    genes_list %<>% filter(genes %in% rownames(seuratObj@data))
  }

  if (translate_gene_names) {
    seuratObj <- correctGeneNames(seuratObj)
    genes_plot <- checkGeneSymbols(
      x = genes_plot,
      unmapped.as.na = FALSE
    ) %>%
      dplyr::pull(Suggested.Symbol) %>%
      unique()
  }

  genes_not_found <- genes_plot %>%
    as_tibble() %>%
    dplyr::filter(!value %in% rownames(seuratObj@data)) %>%
    dplyr::pull(value) %>%
    unique()
  print(glue("The following genes were not found: {genes_not_found}"))
  genes_plot <- genes_plot %>%
    as_tibble() %>%
    dplyr::filter(value %in% rownames(seuratObj@data)) %>%
    dplyr::pull(value) %>%
    unique()

  ident <- as.factor(x = seuratObj@ident)
  if (group_by != "ident") {
    ident <- as.factor(x = FetchData(
      object = seuratObj,
      vars.all = group_by
    )[, 1])
  }

  data_to_plot <- FetchData(
    object = seuratObj,
    vars.all = genes_plot,
    use.scaled = use_scaled
  ) %>%
    as.data.frame()
  data_to_plot$ident <- ident
  data_to_plot <- rownames_to_column(df = data_to_plot, var = "cell")

  data_to_plot %<>% gather(
    key = genes_plot,
    value = expression, -c(cell, ident)
  )

  data_to_plot %<>%
    group_by(ident, genes_plot) %>%
    summarize(
      avg_exp = mean(expm1(x = expression)),
      pct_exp = PercentAbove(x = expression, threshold = 0),
      n = n()
    )

  data_to_plot$genes_plot <- sub(
    x = data_to_plot$genes_plot,
    pattern = "\\.",
    replacement = "-"
  )
  avg_expr <- AverageExpression(
    object = SetAllIdent(seuratObj, group_by),
    genes.use = genes_plot,
    show.progress = FALSE
  ) %>%
    scale()

  if (!is.null(filter_exp_pct)) {
    avg_detect <- AverageDetectionRate(
      object = seuratObj,
      thresh.min = filter_exp_pct_thresh
    )
    avg_detect$highest <- avg_detect %>%
      apply(., MARGIN = 1, FUN = max)
    avg_detect %<>%
      rownames_to_column("gene_name") %>%
      filter(highest > filter_exp_pct) %>%
      dplyr::select(gene_name)
    data_to_plot %<>% filter(genes_plot %in% avg_detect$gene_name)
  }

  if (isTRUE(cluster_x)) {
    gene_dendro <- avg_expr %>%
      dist() %>%
      hclust() %>%
      as.dendrogram()

    data_to_plot$genes_plot <- factor(data_to_plot$genes_plot,
      levels = labels(gene_dendro),
      ordered = TRUE
    )
  }

  if (isTRUE(cluster_y)) {
    id_dendro <- avg_expr %>%
      t() %>%
      dist() %>%
      hclust() %>%
      as.dendrogram()

    data_to_plot$ident <- factor(data_to_plot$ident,
      levels = labels(id_dendro),
      ordered = TRUE
    )
  }

  data_to_plot %<>%
    ungroup() %>%
    group_by(genes_plot) %>%
    mutate(avg_exp_scale = compositions::normalize(x = avg_exp))

  data_to_plot %<>%
    group_by(genes_plot) %>%
    filter(max(avg_exp_scale) > filter_exp_level)

  if (!isTRUE(cluster_y)) {
    data_to_plot <- data_to_plot[mixedorder(data_to_plot$ident), ]
  }
  if (!isTRUE(cluster_x) & isTRUE(order_genes)) {
    data_to_plot <- data_to_plot[mixedorder(data_to_plot$genes_plot), ]
  } else {
    data_to_plot <- data_to_plot[genes_plot, ]
  }

  if (annotated_gene_list) {
    data_to_plot$annotations <- plyr::mapvalues(
      x = data_to_plot$genes_plot,
      from = genes_list$genes,
      to = as.character(genes_list$annotations)
    )
  }

  g <- data_to_plot %>%
    ggplot(aes(
      x = genes_plot,
      y = ident,
      size = pct_exp,
      color = avg_exp_scale
    )) +
    geom_point() +
    theme(
      axis.text.x = element_text(angle = x_lab_rot_angle,
                                 hjust = 1,
                                 size = x_lab_size),
      axis.text.y = element_text(size = y_lab_size)
    ) +
    labs(x = x_axis_title,
         y = y_axis_title,
         size = pct_legend_title,
         color = scale_legend_title) +
    scale_radius(range = c(0, 5))

  if (!is.null(colors_use)) {
    g <- g + scale_color_gradientn(colors = make_color_scale(
      palette = colors_use,
      gradations = 100
    ))
  } else {
    g <- g + scale_color_continuous(low = "#EEEEEE", high = "#FF0000")
  }

  if (annotated_gene_list) {
    g <- g + facet_grid(
      cols = vars(annotations),
      scales = "free_x",
      space = "free_x"
    ) +
      theme(strip.text.x = element_text(size = x_lab_size))
  }

  if (isTRUE(do_return)) {
    return(g)
  } else {
    g
  }
}
