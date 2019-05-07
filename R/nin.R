#' Negative value matching 
#'
#' @name %nin%
#' @rdname negative_match
#' @keywords internal
#' @export
#' @importFrom purrr compose
`%nin%` <- compose(`!`, `%in%`)