#' GObubbleplot
#'
#' Produces a Bubble Plot for the genes of a given GO term.
#'
#' @param object Seurat object
#' @param go_term Gene Ontology term identifier (i.e. GO:0009615)
#' @param group_by Factor by which to group cells.  (default: ident)
#' @param gene_filter A list of gene names to filter the GO term members against.
#' (default: all genes in object)
#' @param filter_exp_pct Display only genes that are expressed above this
#' fraction of cells in at least one group. (default: NULL)
#' @param filter_exp_pct_thresh Threshold for expression fraction. (default: 0)
#' @param filter_exp_level Display only genes that are expressed above this
#' level in at least one group. (default: 0)
#' @param x_lab_size x-axis label font size. Default: 9
#' @param y_lab_size y-axis label font size. Default: 9
#' @param x_lab_rot_angle x-axis label rotation. Default: 45
#' @param cluster_x Should x-axis variables be clustered? Default: TRUE
#' @param cluster_y Should x-axis variables be clustered? Default: TRUE
#' @param colors_use Colors to use for coloring gene expression
#' @param ... Additional options to pass to bubbleplot
#'
#' @import org.Hs.eg.db
#' @importFrom BiocGenerics get mget
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
GObubbleplot <-
  function(
    object,
    go_term,
    group_by              = "ident",
    gene_filter           = NULL,
    filter_exp_pct        = NULL,
    filter_exp_pct_thresh = 0,
    filter_exp_level      = 0,
    x_lab_size            = 9,
    y_lab_size            = 9,
    x_lab_rot_angle       = 45,
    cluster_x             = TRUE,
    cluster_y             = TRUE,
    colors_use            = NULL,
    ...
    ){

  if (is.null(gene_filter)) {
    gene_filter <- rownames(object)
  }

  go_genes_to_plot <- unlist(mget(get(go_term,
                                      org.Hs.egGO2ALLEGS),
                                  org.Hs.egSYMBOL))
  go_genes_to_plot <- go_genes_to_plot[which(go_genes_to_plot %in% gene_filter)]

  if (length(go_genes_to_plot) > 0) {
    gg <- bubbleplot(
      object,
      features_plot         = unique(go_genes_to_plot),
      group_by              = group_by,
      filter_exp_pct        = filter_exp_pct,
      filter_exp_pct_thresh = filter_exp_pct_thresh,
      filter_exp_level      = filter_exp_level,
      x_lab_size            = x_lab_size,
      y_lab_size            = y_lab_size,
      x_lab_rot_angle       = x_lab_rot_angle,
      cluster_x             = cluster_x,
      cluster_y             = cluster_y,
      colors_use            = colors_use,
      ...)
  } else {
    stop("No genes for that term are expressed in the dataset.")
  }

  gg
}
