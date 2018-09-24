#' correctGeneName
#'
#' Convert the gene names in a Seurat object to the accepted HUGO gene names
#'
#' @param seuratObj
#'
#' @importFrom HGNChelper checkGeneSymbols
#' @importFrom plyr mapvalues
#'
#' @return seurat object
#' @export
#'
#' @examples
correctGeneNames <- function(seuratObj){
  corrected.names <- checkGeneSymbols(rownames(seuratObj@raw.data),
                                      unmapped.as.na = FALSE)

  rownames(seuratObj@raw.data) <- mapvalues(x = rownames(seuratObj@raw.data),
                                            from = corrected.names$x,
                                            to = corrected.names$Suggested.Symbol)

  rownames(seuratObj@data) <- mapvalues(x = rownames(seuratObj@data),
                                        from = corrected.names$x,
                                        to = corrected.names$Suggested.Symbol)

  rownames(seuratObj@scale.data) <- mapvalues(x = rownames(seuratObj@scale.data),
                                              from = corrected.names$x,
                                              to = corrected.names$Suggested.Symbol)

  return(seuratObj)
}
