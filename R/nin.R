#' Negative value matching 
#'
#' @name %nin%
#' @rdname negative_match
#' @keywords internal
#' @importFrom purrr compose
`%nin%` <- compose(`!`, `%in%`)