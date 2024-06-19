#' Appends distance between interaction anchors 
#' 
#' Appends distance between interaction anchors, using 
#' `InteractionSet::pairdist`
#'
#' @param x The query GInteractions
#' @param type A character string specifying the type of distance to compute. Can take values of "mid", "gap", "span", "diag" or "intra".
#' @param colname name of column to hold pair distance values
#' 
#' @return The GInteractions with an additional column containing the 
#'   distance between each pair of anchors.
#' 
#' @rdname add-pairdist
#' 
#' @export
#'
#' @examples
#' gi <- read.table(text = "
#' chr1 100 200 chr1 5000 5100 bedpe_example1 30 + -
#' chr1 1000 5000 chr2 3000 3800 bedpe_example2 100 + -",
#' col.names = c(
#'   "seqnames1", "start1", "end1", 
#'   "seqnames2", "start2", "end2", "name", "score", "strand1", "strand2")
#' ) |> as_ginteractions()
#' 
#' add_pairdist(gi)
#' @export

add_pairdist <- function(x, type = 'mid', colname = 'pairdist') {

    if (colname %in% names(GenomicRanges::mcols(x))){
        stop(paste0(colname, " already exists in destination metadata"))
    }
  
    if (is.null(GenomicRanges::mcols(x))){
        # handle IRanges NULL adding X column of NA's
        meta <- S4Vectors::DataFrame("distance" = NA_integer_)
        names(meta) <- colname
        GenomicRanges::mcols(x) <- meta
    } else {
        GenomicRanges::mcols(x)[[colname]] <-  NA_integer_
    }

    GenomicRanges::mcols(x)[[colname]] <- InteractionSet::pairdist(x, type)

    x
}
