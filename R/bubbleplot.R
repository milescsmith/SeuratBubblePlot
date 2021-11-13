#' @title bubbleplot
#' 
#' @description Display percentage of cells expressing and levels for a set of features.
#'
#' @param object Seurat object
#' @param features_plot Either a list of features or a data frame of annotated features
#' to display (see 'annotated_feature_list'). Note: Gene and protein names
#' may be converted to the proper feature name automagically by
#' HGNChelper::checkGeneSymbols if 'translate_feature_names' is TRUE.
#' Features not appearing in the dataset are skipped.
#' @param assay Assay to plot. Default: "RNA"
#' @param slot Slot to plot. Default: "data"
#' @param grouping_var Variable by which to group cells.  Can be cluster identities
#' or any column from the meta.data slot. Default: ident
#' @param filter_exp_pct Display only features that are expressed above this
#' fraction.  The method by which and groups in which this is judged are controlled by
#' `filter_apply_method` and `filter_apply_group`, respectively.  Default: NULL
#' @param filter_apply_method How should `filter_exp_pct` be applied - should `any` or `all` groups
#' be above the desired expression level? Default: "any"
#' @param filter_apply_group What identity groups should be considered when applying 
#' `filter_exp_pct?  Only works with the active identity or whatever is specified in `grouping_var`. 
#' Default: NULL
#' @param filter_exp_pct_thresh Threshold for expression fraction. Default: 0
#' @param filter_exp_level Display only features that are expressed above this
#' level in at least one group. Default: 0
#' @param avg_func Base function to call for measuring average expression. Default: "mean"
#' @param x_lab_size Font size for the x-axis labels. Default: 9
#' @param y_lab_size Font size for the y-axis labels. Default: 9
#' @param x_lab_rot_angle Angle to rotate the x-axis labels. Default: 45°
#' @param preserve_feature_order Should the features by displayed in order in which
#' they are given? Overrides cluster_x. Default: FALSE.
#' @param cluster_x Arrange the x-axis variables using hierarchical clustering.
#' Default: TRUE
#' @param cluster_y Arrange the y-axis variables using hierarchical clustering.
#' Default: FALSE
#' @param colors_use Color palette to use to display expression levels.
#' Default: "Reds"
#' @param translate_feature_names Should feature names be checked and translated to
#' their correct HUGO Gene Nomenclature Committee name? Default: FALSE
#' @param annotated_feature_list Is the feature list annotated?  If so, the features will be
#' faceted using their annontations.  Requires that 'features_plot' is a two column
#' with the annotations in a column named 'annotations' and features in a column
#' named 'features'. Default: FALSE
#' @param do_return Return a ggplot2 object instead of displaying
#' @param verbose Show extra output information, like the features that were not 
#' found? Default: FALSE
#' @param x_axis_title x-axis title
#' @param y_axis_title y-axis title
#' @param pct_legend_title Percent expressed legend title
#' @param scale_legend_title Scale legend title.
#' @param ... ignored
#'
#' @import ggplot2
#' @import Seurat
#' @importFrom dplyr group_by summarise mutate ungroup select pull filter recode 
#' summarise_if mutate_at top_n n
#' @importFrom tibble rownames_to_column as_tibble column_to_rownames
#' @importFrom tidyr gather pivot_longer
#' @importFrom stats hclust dist as.dendrogram order.dendrogram
#' @importFrom compositions normalize
#' @importFrom gtools mixedorder
#' @importFrom HGNChelper checkGeneSymbols
#' @importFrom glue glue
#'
#' @return if isTRUE(do_return), a ggplot2 object
#' @export
#'
#' @examples
bubbleplot <- function(object, ...){
  UseMethod("bubbleplot")
}

#' @rdname bubbleplot
#' @method bubbleplot Seurat
#' @return
#' @export
bubbleplot.Seurat <- function(object,
                              features_plot,
                              assay = NULL,
                              slot = "data",
                              filter_exp_pct = NULL,
                              filter_apply_method = "any",
                              filter_apply_group = NULL,
                              filter_exp_pct_thresh = 0,
                              filter_exp_level = 0,
                              avg_func = "mean",
                              grouping_var = "ident",
                              x_lab_size = 9,
                              y_lab_size = 9,
                              x_axis_title = "Features",
                              y_axis_title = "Grouping",
                              pct_legend_title = "Percent group expressing",
                              scale_legend_title = "Average scaled expression",
                              x_lab_rot_angle = 45,
                              preserve_feature_order = FALSE,
                              cluster_x = TRUE,
                              cluster_y = FALSE,
                              colors_use = NULL,
                              translate_feature_names = FALSE,
                              annotated_feature_list = FALSE,
                              do_return = FALSE,
                              verbose = FALSE,
                              ...) {
  if (isTRUE(annotated_feature_list)) {
    features_list <- features_plot
    features_plot <- features_list$features
    if (is.null(assay)) {
      features_list %<>% filter(features %in% rownames(object))
    }else{
      features_list %<>% filter(features %in% rownames(object[[assay]]))
    }
  }

  if (isTRUE(translate_feature_names)) {
    object <- correctGeneNames(object)
    features_plot <- checkGeneSymbols(x = features_plot,
                                      unmapped.as.na = FALSE) %>%
      pull(Suggested.Symbol) %>%
      unique()
  }
  original_features_order <- features_plot

  ident <- as.factor(x = Idents(object))
  if (grouping_var != "ident") {
    Idents(object) <- grouping_var
    ident <- as.factor(x = FetchData(object = object,
                                     vars = grouping_var,
                                     slot = slot)[, 1])
  }

  if (!is.null(assay)){
    features_plot <- glue("{tolower(assay)}_{features_plot}") %>% as.character()
  }
  data_to_plot <- FetchData(object = object,
                            vars = features_plot) %>%
    as_tibble(rownames = "cell")
  data_to_plot$ident <- ident

  data_to_plot %<>% pivot_longer(cols = -c(cell, ident),
                                 names_to = "features_plot",
                                 values_to = "expression")
  
  features_not_found <- features_plot[features_plot %nin% data_to_plot[["features_plot"]]] %>% unique()
  features_plot <- features_plot[features_plot %in% data_to_plot[["features_plot"]]] %>% unique()
  
  if (isTRUE((verbose))){
    message(glue("The following features were not found: {features_not_found}")) 
  }

  data_to_plot %<>%
    group_by(ident, features_plot) %>%
    summarise(avg_exp = do.call(what = avg_func, args = list(expm1(x = expression))),
              pct_exp = PercentAbove(x = expression, threshold = filter_exp_pct_thresh),
              n = n())

  avg_expr <- FetchData(object = object, 
                       vars = c(features_plot, grouping_var)) %>% 
    as_tibble(rownames = "cell") %>% 
    group_by(ident = get(grouping_var)) %>% 
    summarise_if(is.numeric, get(x = avg_func)) %>% 
    mutate_at(features_plot, normalize) %>%
    as.data.frame() %>% 
    column_to_rownames("ident") %>%
    t()
  
  if (!is.null(filter_exp_pct)) {
    avg_detect <- DetectionRate(object = object,
                                assay = assay,
                                features = features_plot,
                                thresh.min = filter_exp_pct_thresh) %>% 
      as_tibble(rownames="features") %>% 
      pivot_longer(-features, names_to = "ident")
    
    if (!is.null(filter_apply_group)){
      avg_detect %<>% 
        filter(ident %in% filter_apply_group)
    }
    
    if (filter_apply_method == "all"){
      thresholded_features <- avg_detect %>% 
        group_by(features) %>% 
        summarise(minimum = min(value)) %>% 
        filter(minimum >= filter_exp_pct) %>%
        pull(features)
    } else if (filter_apply_method == "any") {
      thresholded_features <- avg_detect %>% 
        group_by(features) %>% 
        top_n(1, value) %>% 
        filter(value >= filter_exp_pct)%>%
        pull(features)
    } else {
      stop("That is not a valid filtering method")
    }
    
    data_to_plot %<>% filter(features_plot %in% thresholded_features)
  }

  if (isTRUE(cluster_x)) {
    feature_dendro <- avg_expr %>%
      dist() %>%
      hclust() %>%
      as.dendrogram()

    data_to_plot$features_plot <- factor(data_to_plot$features_plot,
                                         levels = labels(feature_dendro),
                                         ordered = TRUE)
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
    group_by(features_plot) %>%
    mutate(avg_exp_scale = normalize(x = avg_exp))

  data_to_plot %<>%
    group_by(features_plot) %>%
    filter(max(avg_exp_scale) > filter_exp_level)

  if (!isTRUE(cluster_y)) {
    data_to_plot <- data_to_plot[mixedorder(data_to_plot$ident), ]
  }
  if (!isTRUE(cluster_x)) {
    data_to_plot <- data_to_plot[mixedorder(data_to_plot$features_plot), ]
  }
  
  if (!is.null(assay)){
    data_to_plot$features_plot %<>% str_remove(pattern = glue("{tolower(assay)}_"))
  }
  if (isTRUE(preserve_feature_order)){
    data_to_plot$features_plot <- factor(data_to_plot$features_plot,
                                         levels = unique(original_features_order))
  }
  
  if (annotated_feature_list) {
    rename_list <- as.character(features_list$annotations)
    names(rename_list) <- features_list$features
    data_to_plot$annotations <- recode(.x = data_to_plot$features_plot,
                                       !!!rename_list)
  }
  
  g <- data_to_plot %>%
    ggplot(aes(x = features_plot,
               y = ident,
               size = pct_exp,
               color = avg_exp_scale)) +
    geom_point() +
    theme(axis.text.x = element_text(angle = x_lab_rot_angle,
                                     hjust = 1,
                                     size = x_lab_size),
          axis.text.y = element_text(size = y_lab_size)) +
    labs(x = x_axis_title,
         y = y_axis_title,
         size = pct_legend_title,
         color = scale_legend_title) +
    scale_radius(range = c(0, 5))

  if (!is.null(colors_use)) {
    g <- g + scale_color_gradientn(colors = make_color_scale(palette = colors_use,
                                                             gradations = 100),
                                   limits = c(0,1))
  } else {
    g <- g + scale_color_continuous(low = "#EEEEEE",
                                    high = "#FF0000",
                                    limits = c(0,1))
  }

  if (annotated_feature_list) {
    g <- g + facet_grid(cols = vars(annotations),
                        scales = "free_x",
                        space = "free_x") +
      theme(strip.text.x = element_text(size = x_lab_size))
  }

  if (isTRUE(do_return)) {
    return(g)
  } else {
    g
  }
}
