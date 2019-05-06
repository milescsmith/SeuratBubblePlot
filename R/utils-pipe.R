#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' Compound Assignment Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%<>\%}} for details.
#'
#' @name %<>%
#' @rdname compound_pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %<>%
#' @usage lhs \%<>\% rhs
NULL

#' Negative value matching 
#'
#' @name %nin%
#' @rdname negative_match
#' @keywords internal
#' @export
#' @importFrom purrr compose
#' @usage lhs \%nin\% rhs
`%nin%` <- compose(`!`, `%in%`)