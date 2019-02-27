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
correctGeneNames <- function(object, ...) {
  UseMethod("correctGeneNames")
}

#' @rdname correctGeneNames
#' @method correctGeneNames seurat
#' @export
#' @return
correctGeneNames.seurat <- function(object) {
  corrected_names <- checkGeneSymbols(rownames(object@raw.data),
    unmapped.as.na = TRUE) %>%
    filter(is.na(Suggested.Symbol))

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
    unmapped.as.na = TRUE) %>%
    filter(is.na(Suggested.Symbol))

  corrected_names[["Suggested.Symbol"]] %<>% make.names(unique = TRUE)
  sn <- c(
    "counts",
    "data",
    "scale.data",
    "meta.features"
  )
  for (i in sn) {
    suppressMessages(
      rownames(slot(object@assays[[assay]], i)) <- mapvalues(
        x = rownames(slot(object@assays[[assay]], i)),
        from = corrected_names[["x"]],
        to = as.character(corrected_names[["Suggested.Symbol"]])
      )
    )
  }
  suppressMessages(
    object@assays[[assay]]@var.features <- mapvalues(
      x = object@assays[[assay]]@var.features,
      from = corrected_names[["x"]],
      to = as.character(corrected_names[["Suggested.Symbol"]])
    )
  )

  return(object)
}
