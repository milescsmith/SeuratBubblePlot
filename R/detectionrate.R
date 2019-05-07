#' @title DetectionRate
#' @Description Probability of detection by identity class
#' Stolen from Seurat 2.3.4 (and rewritten), since it was removed from Seurat 3
#'
#' For each gene, calculates the probability of detection for each identity
#' class.
#'
#' @param object Seurat object
#' @param thresh.min Minimum threshold to define 'detected' (log-scale)
#' @param features Which features to calculate detection rate for. Default: NULL (= all)
#' @param slot_use Slot to pull data from.  Default: "data"
#' @param ... ignored
#'
#' @importFrom stringr str_remove
#' @importFrom purrr map_dfr map
#' @importFrom glue glue
#'
#' @return Returns a matrix with genes as rows, identity classes as columns.
#'
#' @export
#'
#' @examples
DetectionRate <- function(object, ...){
  UseMethod("DetectionRate")
}

#' @rdname DetectionRate
#' @method DetectionRate Seurat
#' @importFrom Seurat FetchData WhichCells Idents
#' @export
#' @return
DetectionRate.Seurat <- function(object,
                                 features = NULL,
                                 slot_use = "data",
                                 thresh.min = 0,
                                 ...) {
  assay <- DefaultAssay(object)
  ident_use <- Idents(object)
  data_all <- map_dfr(sort(x = unique(x = ident_use)), 
                      function(i) {
                        temp_cells <- WhichCells(object = object, 
                                                 ident = i)
                        vars_use <- glue("{tolower(assay)}_{features}") %>% 
                          as.character()
                        data.temp <- map(FetchData(object, 
                                                   vars = vars_use, 
                                                   cells = temp_cells, 
                                                   slot = slot_use), 
                                         function(x){
                                           sum(x > thresh.min)/length(x = x)
                                           }) 
                        }) %>% 
    t()
  colnames(x = data_all) <- sort(x = unique(x = ident_use))
  rownames(x = data_all) %<>% str_remove(glue("{tolower(assay)}_"))
  return(data_all)
}