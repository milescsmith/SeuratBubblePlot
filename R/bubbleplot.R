#' bubbleplot
#'
#' Display percentage of cells expressing and levels for a set of genes.
#'
#' @param seuratObj Seurat object
#' @param genes.plot Either a list of genes or a data frame of annotated genes
#'   to display (see 'annotated.gene.list'). Note: Gene and protein names
#'   may be converted to the proper gene name automagically by
#'   HGNChelper::checkGeneSymbols if 'translate.gene.names' is TRUE.
#'   Genes not appearing in the dataset are skipped.
#' @param filter.exp.pct Display only genes that are expressed above this
#'   fraction of cells in at least one group. (default: NULL)
#' @param filter.exp.pct.thresh Threshold for expression fraction. (default: 0)
#' @param filter.exp.level Display only genes that are expressed above this
#'   level in at least one group. (default: 0)
#' @param group.by Variable by which to group cells.  Can be cluster identities
#'   or any column from the meta.data slot. (default: ident)
#' @param x.lab.size Font size for the x-axis labels. (default: 9)
#' @param y.lab.size Font size for the y-axis labels. (default: 9)
#' @param x.lab.rot.angle Angle to rotate the x-axis labels. (default: 45°)
#' @param cluster.x Arrange the x-axis variables using hierarchical clustering.
#'   (default: TRUE)
#' @param cluster.y Arrange the y-axis variables using hierarchical clustering.
#'   (default: FALSE)
#' @param colors.use Color palette to use to display expression levels.
#'   (default: "Reds")
#' @param translate.gene.names Should gene names be checked and translated to
#'   their correct HUGO Gene Nomenclature Committee name? (default: FALSE)
#' @param annotated.gene.list Is the gene list annotated?  If so, the genes will be
#'   faceted using their annontations.  Requires that 'genes.plot' is a two column
#'   with the annotations in a column named 'annotations' and genes in a column
#'   named 'genes'. (default: FALSE)
#' @param do.return Return a ggplot2 object instead of displaying
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
#' @return if isTRUE(do.return), a ggplot2 object
#' @export
#'
#' @examples BubblePlot(seuratObj = obj, genes.plot = c("IFIT1","IFITM1","IFITM3"), group.by = "treatment")
bubbleplot <- function(seuratObj,
                       genes.plot,
                       use.scaled = FALSE,
                       filter.exp.pct = NULL,
                       filter.exp.pct.thresh = 0,
                       filter.exp.level = 0,
                       group.by = "ident",
                       x.lab.size = 9,
                       y.lab.size = 9,
                       x.axis.title = "Genes",
                       y.axis.title = "Grouping",
                       pct.legend.title = "Percent group expressing",
                       scale.legend.title = "Average scaled expression",
                       x.lab.rot.angle = 45,
                       cluster.x = TRUE,
                       cluster.y = FALSE,
                       colors.use = NULL,
                       translate.gene.names = FALSE,
                       annotated.gene.list = FALSE,
                       do.return = FALSE) {

  if (annotated.gene.list){
    genes.list <- genes.plot
    genes.plot <- genes.list$genes
    genes.list %<>% filter(genes %in% rownames(seuratObj@data))
  }

  if (translate.gene.names) {
    genes.plot <- checkGeneSymbols(
      x = genes.plot,
      unmapped.as.na = FALSE
    ) %>%
      dplyr::pull(Suggested.Symbol) %>%
      unique()
  }

  genes.not.found <- genes.plot %>%
    as_tibble() %>%
    dplyr::filter(!value %in% rownames(seuratObj@data)) %>%
    dplyr::pull(value) %>%
    unique()
  print(glue("The following genes were not found: {genes.not.found}"))
  genes.plot <- genes.plot %>%
    as_tibble() %>%
    dplyr::filter(value %in% rownames(seuratObj@data)) %>%
    dplyr::pull(value) %>%
    unique()

  ident <- as.factor(x = seuratObj@ident)
  if (group.by != "ident") {
    ident <- as.factor(x = FetchData(
      object = seuratObj,
      vars.all = group.by
    )[, 1])
  }

  seuratObj <- correctGeneNames(seuratObj)
  data.to.plot <- FetchData(
    object = seuratObj,
    vars.all = genes.plot,
    use.scaled = use.scaled) %>%
    as.data.frame()
  data.to.plot$ident <- ident
  data.to.plot <- rownames_to_column(df = data.to.plot, var = "cell")

  data.to.plot <- data.to.plot %>% gather(
    key = genes.plot,
    value = expression, -c(cell, ident)
  )

  data.to.plot <- data.to.plot %>%
    group_by(ident, genes.plot) %>%
    summarize(
      avg.exp = mean(expm1(x = expression)),
      pct.exp = PercentAbove(x = expression, threshold = 0),
      n = n()
    )

  data.to.plot$genes.plot <- sub(
    x = data.to.plot$genes.plot,
    pattern = "\\.",
    replacement = "-"
  )
  avg.expr <- AverageExpression(
    object = SetAllIdent(seuratObj, group.by),
    genes.use = genes.plot,
    show.progress = FALSE
  ) %>%
    scale()

  if (!is.null(filter.exp.pct)) {
    avg.detect <- AverageDetectionRate(
      object = seuratObj,
      thresh.min = filter.exp.pct.thresh
    )
    avg.detect$highest <- avg.detect %>%
      apply(., MARGIN = 1, FUN = max)
    avg.detect %<>%
      rownames_to_column("gene_name") %>%
      filter(highest > filter.exp.pct) %>%
      dplyr::select(gene_name)
    data.to.plot %<>% filter(genes.plot %in% avg.detect$gene_name)
  }

  if (isTRUE(cluster.x)) {
    gene.dendro <- avg.expr %>%
      dist() %>%
      hclust() %>%
      as.dendrogram()

    data.to.plot$genes.plot <- factor(data.to.plot$genes.plot,
                                      levels = labels(gene.dendro),
                                      ordered = TRUE
    )
  }

  if (isTRUE(cluster.y)) {
    id_dendro <- avg.expr %>%
      t() %>%
      dist() %>%
      hclust() %>%
      as.dendrogram()

    data.to.plot$ident <- factor(data.to.plot$ident,
                                 levels = labels(id_dendro),
                                 ordered = TRUE
    )
  }

  data.to.plot <- data.to.plot %>%
    ungroup() %>%
    group_by(genes.plot) %>%
    mutate(avg.exp.scale = compositions::normalize(x = avg.exp))

  data.to.plot %<>%
    group_by(genes.plot) %>%
    filter(max(avg.exp.scale) > filter.exp.level)

  if (!isTRUE(cluster.y)) {
    data.to.plot <- data.to.plot[mixedorder(data.to.plot$ident), ]
  }
  if (!isTRUE(cluster.x)) {
    data.to.plot <- data.to.plot[mixedorder(data.to.plot$genes.plot), ]
  }

  if(annotated.gene.list) {
    data.to.plot$annotations <- plyr::mapvalues(x = data.to.plot$genes.plot,
                                               from = genes.list$genes,
                                               to = as.character(genes.list$annotations))
  }

  g <- data.to.plot %>%
    ggplot(aes(
      x = genes.plot,
      y = ident,
      size = pct.exp,
      color = avg.exp.scale
    )) +
    geom_point() +
    theme(
      axis.text.x = element_text(angle = x.lab.rot.angle, hjust = 1, size = x.lab.size),
      axis.text.y = element_text(size = y.lab.size)
    ) +
    labs(x = x.axis.title, y = y.axis.title, size = pct.legend.title, color = scale.legend.title) +
    scale_radius(range = c(0, 5))

  if (!is.null(colors.use)) {
    g <- g + scale_color_gradientn(colors = make_color_scale(palette = colors.use,
                                                             gradations = 100))
  } else {
    g <- g + scale_color_continuous(low = "#EEEEEE", high = "#FF0000")
  }

  if(annotated.gene.list) {
    g <- g + facet_grid(cols = vars(annotations),
                        scales = "free_x",
                        space = "free_x") +
      theme(strip.text.x = element_text(size = x.lab.size))
  }

  if (isTRUE(do.return)) {
    return(g)
  } else {
    g
  }
}
