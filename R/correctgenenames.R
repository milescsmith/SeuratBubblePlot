#' @title correctGeneNames
#'
#' @description Convert the gene names in a Seurat object to the accepted HUGO gene names
#'
#' @param object Seurat object
#' @param assay assay to pull data from. If NULL, the default assay is used.  Default: NULL
#' @param ... ignored
#'
#' @import Seurat
#' @importFrom HGNChelper checkGeneSymbols
#' @importFrom dplyr recode
#' @importFrom methods slot<- slot
#'
#' @return seurat object
#' @export
#'
#' @examples
correctGeneNames <- function(object, ...) {
  UseMethod("correctGeneNames")
}


#' @rdname correctGeneNames
#' @method correctGeneNames Seurat
#' @export
#' @return
correctGeneNames.Seurat <-
  function(object,
           assay = NULL,
           eliminate_unknown = FALSE,
           ...) {
    if (is.null(assay)) {
      assay <- DefaultAssay(object)
    }

    corrected_names <- suppressWarnings(HGNChelper::checkGeneSymbols(rownames(object),
      unmapped.as.na = TRUE
    )) |> 
      filter(is.na(Suggested.Symbol))

    corrected_names[["Suggested.Symbol"]] <- corrected_names[["Suggested.Symbol"]] |> make.names(unique = TRUE)
    sn <- c(
      "counts",
      "data",
      "scale.data",
      "meta.features"
    )
    rename_list <- as.character(corrected_names[["Suggested.Symbol"]])
    names(rename_list) <- corrected_names[["x"]]
    for (i in sn) {
      suppressMessages(
        rownames(slot(object@assays[[assay]], i)) <- recode(
          .x = rownames(slot(object@assays[[assay]], i)),
          !!!rename_list
        )
      )
    }
    suppressMessages(
      object@assays[[assay]]@var.features <- recode(
        .x = object@assays[[assay]]@var.features,
        !!!rename_list
      )
    )

    if (isTRUE(eliminate_unknown)) {
      object <- object[
        stringr::str_detect(
          string = rownames(object),
          pattern = "^NA.",
          negate = TRUE
        ),
      ]
    }


    object
  }
