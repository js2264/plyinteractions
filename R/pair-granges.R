#' Pairwise combination of a GRanges object
#' 
#' Create a GInteractions object from a GRanges object, 
#' containing all possible entry pairs
#'
#' @param x A GRanges object
#' 
#' @return A GInteractions object
#' 
#' @rdname pair-granges
#' 
#' @export
#'
#' @examples
#' gr <- read.table(text = "
#' chr1 100 200
#' chr1 5000 5100
#' chr1 1000 5000
#' chr2 3000 3800",
#' col.names = c(
#'   "seqnames", "start", "end" 
#' )) |> plyranges::as_granges()
#' 
#' pair_granges(gr)
#' @export

pair_granges <- function(x) {

    combs <- combn(length(x), 2)
    InteractionSet::GInteractions(combs[1,], combs[2,], gr)

}
