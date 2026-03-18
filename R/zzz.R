#' plyinteractions: Extending tidyomics verbs to genomic interactions
#'
#' @description
#' plyinteractions verbs treat `GInteractions` objects as tabular data using 
#'  `dplyr`-like verbs. The functions and methods in `plyinteractions`
#'  provide a grammatical approach to manipulate `GInteractions`, to 
#'  facilitate their integration in genomic analysis workflows.
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/js2264/plyinteractions}
#'   \item Report bugs at \url{https://github.com/js2264/plyinteractions/issues}
#' }
#'
#' @author Jacques Serizay 
#'
#' @docType package
#' @name plyinteractions-package
#' @aliases plyinteractions
#' @include attach.R
#' @keywords internal
"_PACKAGE"

.onAttach <- function(libname, pkgname) {
    attached <- tidyverse_attach()
}
