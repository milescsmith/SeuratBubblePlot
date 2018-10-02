#' PercentAbove
#'
#' Return the percentage of a list that is above a threshold
#' @param x A list of numeric values.
#' @param threshold A numeric threshold value.
#'
#' @return double
#' @export
#'
#' @examples PercentAbove(x = c(1,2,3), threshold = 2)
PercentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}
