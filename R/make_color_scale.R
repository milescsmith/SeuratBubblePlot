#' make_color_scale
#'
#' Using a list of colors, create a gradient of a discrete number of colors.
#'
#' @param palette One of the palettes in RColorBrewer, viridis, or a list of colors understood by colors()
#' @param gradations The number of colors to return
#'
#' @return A list of hexidecimal numbers.
#' @export
#'
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @importFrom viridis viridis
#' @importFrom grDevices colorRampPalette
#'
#' @examples make_color_scale(palette = 'Reds', n = 100)
make_color_scale <- function(palette = "viridis", gradations = 10){
  if (length(palette == 1)){
    if (palette %in% rownames(brewer.pal.info)){
      pal <- colorRampPalette(brewer.pal(brewer.pal.info[palette,]$maxcolors,palette))(gradations)
    } else if (palette %in% c('viridis','inferno','magma','plasma','cividis')){
      pal <- viridis(n = gradations, option = palette)
    }
  } else {
    pal <- colorRampPalette(palette)(gradations)
  }
  return(pal)
}
