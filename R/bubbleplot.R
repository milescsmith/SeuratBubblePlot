#' bubbleplot
#'
#' Display percentage of cells expressing and levels for a set of genes.
#'
#' @param seuratObj Seurat object
#' @param genes.plot A list of genes to display.
#' @param group.by Variable by which to group cells.  Can be cluster identities or any column from the meta.data slot. (default: ident)
#' @param x.lab.size Font size for the x-axis labels. (default: 9)
#' @param y.lab.size Font size for the y-axis labels. (default: 9)
#' @param x.lab.rot.angle Angle to rotate the x-axis labels. (default: 45°)
#' @param clust.x Arrange the x-axis variables using hierarchical clustering. (default: TRUE)
#' @param clust.y Arrange the y-axis variables using hierarchical clustering. (default: TRUE)
#' @param colors.use Color palette to use to display expression levels. (note: not currently implemented)
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
#'
#' @return if isTRUE(do.return), a ggplot2 object
#' @export
#'
#' @examples BubblePlot(seuratObj = obj, genes.plot = c("IFIT1","IFITM1","IFITM3"), group.by = "treatment")
bubbleplot <- function(seuratObj,
                       genes.plot,
                       group.by = 'ident',
                       x.lab.size = 9,
                       y.lab.size = 9,
                       x.lab.rot.angle = 45,
                       clust.x = TRUE,
                       clust.y = TRUE,
                       colors.use = NULL,
                       do.return = FALSE){

  genes.plot <- (genes.plot %>% as_tibble() %>% dplyr::filter(value %in% rownames(seuratObj@data)))$value

  ident <- as.factor(x = seuratObj@ident)
  if (group.by != "ident") {
    ident <- as.factor(x = FetchData(
      object = seuratObj,
      vars.all = group.by
    )[, 1])
  }

  data.to.plot <- data.frame(FetchData(object = seuratObj, vars.all = genes.plot))
  data.to.plot$ident <- ident
  data.to.plot <- rownames_to_column(df = data.to.plot, var = "cell")
  #colnames(data.to.plot)[dim(data.to.plot)[2]] <- "ident"
  data.to.plot <- data.to.plot %>% gather(key = genes.plot,
                                          value = expression, -c(cell, ident))

  data.to.plot <- data.to.plot %>% group_by(ident, genes.plot) %>%
    summarize(avg.exp = mean(expm1(x = expression)),
              pct.exp = PercentAbove(x = expression, threshold = 0),
              n = n())

  data.to.plot$genes.plot <- sub(x = data.to.plot$genes.plot, pattern = "\\.", replacement = "-")
  avg_expr <- AverageExpression(object = SetAllIdent(seuratObj, group.by), genes.use = genes.plot, show.progress = FALSE) %>%
    scale()

  if(isTRUE(clust.x)){
    gene_dendro <- avg_expr %>% dist() %>% hclust %>% as.dendrogram()
    gene_order <- gene_dendro %>% order.dendrogram()
    data.to.plot$genes.plot <- factor(data.to.plot$genes.plot, levels = labels(gene_dendro), ordered = TRUE)
  }

  if(isTRUE(clust.y)){
    id_dendro <- avg_expr %>% t() %>% dist() %>% hclust %>% as.dendrogram()
    id_order <- id_dendro %>% order.dendrogram()
    data.to.plot$ident <- factor(data.to.plot$ident, levels = labels(id_dendro), ordered = TRUE)
  }

  data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>%
    mutate(avg.exp.scale = compositions::normalize(x = avg.exp))

  g <- data.to.plot %>%
    ggplot(aes(x = genes.plot,
               y = ident,
               size = pct.exp,
               color = avg.exp.scale)
           ) +
    geom_point() +
    theme(axis.text.x = element_text(angle=x.lab.rot.angle, hjust = 1, size = x.lab.size),
          axis.text.y = element_text(size = y.lab.size)) +
    scale_radius(range = c(0,5))

  if(!is.null(colors.use)){
    g <- g + scale_color_gradientn(colors = make_color_scale(palette = colors.use, gradations = 100))
  } else {
    g <- g + scale_color_continuous(low = "#EEEEEE",high = "#FF0000")
  }

  if(isTRUE(do.return)){
    return(g)
  } else {
    g
  }
}

#' PercentAbove
#'
#' Return the percentage of a list that is above a threshold
#' @param x A list of numeric values.
#' @param threshold A numeric threshold value.
#'
#' @return double
#' @export
#'
#' @examples PercentAbove(x = c(1,2,3), threshold = 2)
PercentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}


#' GObubbleplot
#'
#' Produces a Bubble Plot for the genes of a given GO term.
#'
#' @param seuratObj Seurat object
#' @param go_term Gene Ontology term identifier (i.e. GO:0046774)
#' @param group.by Factor by which to group cells.  (default: ident)
#' @param filter A list of gene names to filter the GO term members against. (default: all genes in seuratObj)
#' @param do.return If TRUE, return a ggplot2 object instead of displaying chart
#'
#' @import dplyr
#' @import org.Hs.eg.db
#' @importFrom magrittr "%>%"
#'
#' @return if isTRUE(do.return), a ggplot2 object
#' @export
#'
#' @examples GObubbleplot(seuratObj = dataset, go_term = "GO:0002253", filter = dataset@var.genes)
GObubbleplot <- function(seuratObj,
                         go_term,
                         group.by = "ident",
                         filter = NULL,
                         x.lab.size = 9,
                         y.lab.size = 9, 
                         x.lab.rot.angle = 45,
                         clust.x = TRUE, 
                         clust.y = TRUE,
                         colors.use = NULL, 
                         do.return = FALSE){
  if(is.null(filter)){
    filter <- rownames(seuratObj@data)
  }

  go_genes_to_plot <- unlist(mget(get(go_term,org.Hs.egGO2ALLEGS), org.Hs.egSYMBOL))
  go_genes_to_plot <- go_genes_to_plot[which(go_genes_to_plot %in% filter)]

  if(length(go_genes_to_plot) > 0){
    gg <- bubbleplot(seuratObj,
                     genes.plot = unique(go_genes_to_plot),
                     group.by = group.by,
                     x.lab.size = x.lab.size,
                     y.lab.size = y.lab.size, 
                     x.lab.rot.angle = x.lab.rot.angle, 
                     clust.x = clust.x, 
                     clust.y = clust.y, 
                     colors.use = colors.use, 
                     do.return = do.return)
  } else {
    print("No genes for that term are expressed in the dataset.")
    exit()
  }

  if(isTRUE(do.return)){
    return(gg)
  } else {
    gg
  }

}
