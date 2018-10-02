#' GObubbleplot
#'
#' Produces a Bubble Plot for the genes of a given GO term.
#'
#' @param seuratObj Seurat object
#' @param go_term Gene Ontology term identifier (i.e. GO:0046774)
#' @param group_by Factor by which to group cells.  (default: ident)
#' @param filter A list of gene names to filter the GO term members against.
#'   (default: all genes in seuratObj)
#' @param filter_exp_pct Display only genes that are expressed above this
#'   fraction of cells in at least one group. (default: NULL)
#' @param filter_exp_pct_thresh Threshold for expression fraction. (default: 0)
#' @param filter_exp_level Display only genes that are expressed above this
#'   level in at least one group. (default: 0)
#' @param ... Additional options to pass to bubbleplot
#'
#' @param do_return If TRUE, return a ggplot2 object instead of displaying chart
#'
#' @import dplyr
#' @import org.Hs.eg.db
#' @importFrom magrittr "%>%"
#' @importFrom R.utils exit
#'
#' @return if isTRUE(do_return), a ggplot2 object
#' @export
#'
#' @examples GObubbleplot(seuratObj = dataset, go_term = "GO:0002253", filter = dataset@var.genes)
GObubbleplot <- function(seuratObj,
                         go_term,
                         group_by = "ident",
                         filter = NULL,
                         filter_exp_pct = NULL,
                         filter_exp_pct_thresh = 0,
                         filter_exp_level = 0,
                         x_lab_size = 9,
                         y_lab_size = 9,
                         x_lab_rot_angle = 45,
                         cluster_x = TRUE,
                         cluster_y = TRUE,
                         colors_use = NULL,
                         do_return = FALSE,
                         ...){

  if (is.null(filter)) {
    filter <- rownames(seuratObj@data)
  }

  go_genes_to_plot <- unlist(BiocGenerics::mget(
    BiocGenerics::get(go_term,
                      org.Hs.egGO2ALLEGS),
    org.Hs.egSYMBOL))
  go_genes_to_plot <- go_genes_to_plot[which(go_genes_to_plot %in% filter)]

  if (length(go_genes_to_plot) > 0) {
    gg <- bubbleplot(seuratObj,
                     genes_plot = unique(go_genes_to_plot),
                     group_by = group_by,
                     filter.exp_pct = filter_exp_pct,
                     filter_exp_pct_thresh = filter_exp_pct_thresh,
                     filter_exp_level = filter_exp_level,
                     x_lab_size = x_lab_size,
                     y_lab_size = y_lab_size,
                     x_lab_rot_angle = x_lab_rot_angle,
                     cluster_x = cluster_x,
                     cluster_y = cluster_y,
                     colors_use = colors_use,
                     do_return = do_return,
                     ...)
  } else {
    print("No genes for that term are expressed in the dataset.")
    exit()
  }

  if (isTRUE(do_return)) {
    return(gg)
  } else {
    gg
  }

}
