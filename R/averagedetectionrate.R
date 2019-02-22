#' @title DetectionRate
#' @Description Probability of detection by identity class
#' Stolen from Seurat 2.3.4, since it was removed from Seurat 3
#'
#' For each gene, calculates the probability of detection for each identity
#' class.
#'
#' @param object Seurat object
#' @param thresh.min Minimum threshold to define 'detected' (log-scale)
#'
#' @return Returns a matrix with genes as rows, identity classes as columns.
#'
#' @export
#'
#' @examples
#' head(DetectionRate(object = pbmc_small))
DetectionRate <- function(object, ...){
  UseMethod("DetectionRate")
}

#' @rdname DetectionRate
#' @method DetectionRate Seurat
#' @export
#' @return
DetectionRate.Seurat <- function(object, assay = "RNA", thresh.min = 0) {
  DefaultAssay(object) <- assay
  ident.use <- Idents(object)
  data.all <- data.frame(row.names = rownames(object))
  for (i in sort(x = unique(x = ident.use))) {
    temp.cells <- WhichCells(object = object, ident = i)
    data.temp <- apply(
      X = GetAssayData(object, slot = "data", assay = assay)[, temp.cells],
      MARGIN = 1,
      FUN = function(x) {
        return(sum(x > thresh.min)/length(x = x))
      }
    )
    data.all <- cbind(data.all, data.temp)
    colnames(x = data.all)[ncol(x = data.all)] <- i
  }
  colnames(x = data.all) <- sort(x = unique(x = ident.use))
  return(data.all)
}

#' @rdname DetectionRate
#' @method DetectionRate seurat
#' @export
#' @return
DetectionRate.seurat <- function(object, thresh.min = 0) {
  ident.use <- object@ident
  data.all <- data.frame(row.names = rownames(x = object@data))
  for (i in sort(x = unique(x = ident.use))) {
    temp.cells <- WhichCells(object = object, ident = i)
    data.temp <- apply(
      X = object@data[, temp.cells],
      MARGIN = 1,
      FUN = function(x) {
        return(sum(x > thresh.min)/length(x = x))
      }
    )
    data.all <- cbind(data.all, data.temp)
    colnames(x = data.all)[ncol(x = data.all)] <- i
  }
  colnames(x = data.all) <- sort(x = unique(x = ident.use))
  return(data.all)
}
