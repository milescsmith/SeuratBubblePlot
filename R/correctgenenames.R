#' @title correctGeneNames
#'
#' @description Convert the gene names in a Seurat object to the accepted HUGO gene names
#'
#' @param object
#'
#' @import Seurat
#' @importFrom HGNChelper checkGeneSymbols
#' @importFrom plyr mapvalues
#' @importFrom methods slot
#'
#' @return seurat object
#' @export
#'
#' @examples
correctGeneNames <- function(object, ...){
  UseMethod("correctGeneNames")
}

#' @rdname correctGeneNames
#' @method correctGeneNames seurat
#' @export
#' @return
correctGeneNames.seurat <- function(object) {
  corrected_names <- checkGeneSymbols(rownames(object@raw.data),
    unmapped.as.na = FALSE
  )

  rownames(object@raw.data) <- mapvalues(
    x = rownames(object@raw.data),
    from = corrected_names[["x"]],
    to = corrected_names[["Suggested.Symbol"]]
  )

  rownames(object@data) <- mapvalues(
    x = rownames(object@data),
    from = corrected_names[["x"]],
    to = corrected_names[["Suggested.Symbol"]]
  )

  rownames(object@scale.data) <- mapvalues(
    x = rownames(object@scale.data),
    from = corrected_names[["x"]],
    to = corrected_names[["Suggested.Symbol"]]
  )

  return(object)
}


#' @rdname correctGeneNames
#' @method correctGeneNames Seurat
#' @export
#' @return
correctGeneNames.Seurat <- function(object,
                                    assay = "RNA") {
  corrected_names <- checkGeneSymbols(rownames(object),
                                      unmapped.as.na = FALSE
  )

  sn <- c("counts",
          "data",
          "scale.data",
          "meta.features")
  for (i in sn){
    rownames(slot(object@assays[[assay]]),i) <- mapvalues(x = rownames(slot(object@assays[[assay]]),i),
                                                         from = corrected_names[["x"]],
                                                         to = corrected_names[["Suggested.Symbol"]])
  }

  object@data[[assay]]@var.features <- mapvalues(x = object@data[[assay]]@var.features,
                                                         from = corrected_names[["x"]],
                                                         to = corrected_names[["Suggested.Symbol"]])

  return(object)
}
