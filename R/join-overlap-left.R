#' Join overlaps between a query GInteractions and a GRanges
#'
#' @param x A (Pinned)GInteractions object
#' @param y A GRanges object
#' @param maxgap,minoverlap See \code{?\link[GenomicRanges]{countOverlaps}} 
#' in the \pkg{GenomicRanges} package for a description of these arguments
#' @param suffix Suffix to add to metadata 
#' columns (character vector of length 2, default to `c(".x", ".y")`).
#' @param distance If `TRUE`, add a `distance` column: for each interaction 
#' joined to a range of `y`, the distance between this range and the 
#' interaction's anchor (the pinned anchor of a PinnedGInteractions, the 
#' closer of the two anchors otherwise), as computed by 
#' `GenomicRanges::distance()`: `0` when they overlap, positive when they are 
#' within `maxgap` of each other. Interactions joined to no range get `NA`.
#'
#' @importFrom plyranges join_overlap_left
#' @importFrom plyranges join_overlap_left_directed
#' @importFrom InteractionSet countOverlaps
#' 
#' @return A GInteractions object, with the metadata columns of `y` (and 
#' `distance`) added. Interactions overlapping several ranges of `y` are 
#' repeated, and those overlapping none are kept, with missing values.
#'
#' @name ginteractions-join-overlap-left
#' 
#' @examples
#' gi <- read.table(text = "  
#'     chr1 11 20 - chr1 21 30 + 
#'     chr1 11 20 - chr1 51 55 + 
#'     chr1 21 30 - chr1 51 55 + 
#'     chr1 21 30 - chr2 51 60 +",  
#' col.names = c(
#'     "seqnames1", "start1", "end1", "strand1", 
#'     "seqnames2", "start2", "end2", "strand2")
#' ) |> as_ginteractions() |> mutate(id = 1:4, type = 'gi')
#' 
#' gr <- GenomicRanges::GRanges(
#'     c("chr1:20-30:+", "chr2:55-65:-")
#' ) |> mutate(id = 1:2, type = 'gr')
#' 
#' gi
#' 
#' gr
#' 
#' ####################################################################
#' # 1. Join overlaps between GInteractions and a subject GRanges
#' ####################################################################
#' 
#' join_overlap_left(gi, gr)
#' 
#' join_overlap_left_directed(gi, gr)
#' 
#' ####################################################################
#' # 2. Join overlaps between PinnedGInteractions and a subject GRanges
#' ####################################################################
#' 
#' gi |> pin_by("first") |> join_overlap_left(gr)
#' 
#' gi |> pin_by("first") |> join_overlap_left_directed(gr)
#' 
#' gi |> pin_by("second") |> join_overlap_left(gr)
#' 
#' gi |> pin_by("second") |> join_overlap_left_directed(gr)
#' 
#' ####################################################################
#' # 3. Distance between the joined ranges
#' ####################################################################
#' 
#' join_overlap_left(gi, gr, maxgap = 25, distance = TRUE)
#' 
#' gi |> pin_by("second") |> join_overlap_left(gr, maxgap = 25, distance = TRUE)
NULL

#' @rdname ginteractions-join-overlap-left
#' @export
join_overlap_left.PinnedGInteractions <- function(
    x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y"), distance = FALSE
) {

    hits <- findOverlaps(
        y, 
        unpin(x), 
        maxgap = maxgap, 
        minoverlap = minoverlap,
        type = 'any', 
        select = 'all', 
        ignore.strand = TRUE, 
        use.region = switch(
            as.character(pin(x)), 
            "1" = "first",
            "2" = "second"
        )
    )

    .hits_to_left(unpin(x), y, hits, suffix, distance = if (distance) {
        .hits_distance(unpin(x), y, hits, as.integer(pin(x)), TRUE)
    }) 

}

#' @rdname ginteractions-join-overlap-left
#' @export
join_overlap_left.GInteractions <- function(
    x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y"), distance = FALSE
) {

    hits <- findOverlaps(
        y, 
        x,
        maxgap = maxgap, 
        minoverlap = minoverlap,
        type = 'any', 
        select = 'all', 
        ignore.strand = TRUE, 
        use.region = 'both'
    )

    .hits_to_left(x, y, hits, suffix, distance = if (distance) {
        .hits_distance(x, y, hits, c(1L, 2L), TRUE)
    }) 

}

#' @rdname ginteractions-join-overlap-left
#' @export
join_overlap_left_directed.PinnedGInteractions <- function(
    x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y"), distance = FALSE
) {

    hits <- findOverlaps(
        y, 
        unpin(x), 
        maxgap = maxgap, 
        minoverlap = minoverlap,
        type = 'any', 
        select = 'all', 
        ignore.strand = FALSE, 
        use.region = switch(
            as.character(pin(x)), 
            "1" = "first",
            "2" = "second"
        )
    )

    .hits_to_left(unpin(x), y, hits, suffix, distance = if (distance) {
        .hits_distance(unpin(x), y, hits, as.integer(pin(x)), FALSE)
    }) 

}

#' @rdname ginteractions-join-overlap-left
#' @export
join_overlap_left_directed.GInteractions <- function(
    x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y"), distance = FALSE
) {

    hits <- findOverlaps(
        y, 
        x,
        maxgap = maxgap, 
        minoverlap = minoverlap,
        type = 'any', 
        select = 'all', 
        ignore.strand = FALSE, 
        use.region = 'both'
    )

    .hits_to_left(x, y, hits, suffix, distance = if (distance) {
        .hits_distance(x, y, hits, c(1L, 2L), FALSE)
    }) 

}

#' @importFrom S4Vectors decode
#' @importFrom S4Vectors Rle
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
.na_cols <- function(x) {
    if (is(x, "IRanges"))
        return(IRanges::IRanges(0L, -1L))
    if (is(x, "GRanges"))
        return(GenomicRanges::GRanges(".", IRanges::IRanges(0L, -1L)))
    if (is(x, "Rle"))
        x <- S4Vectors::decode(x)
    ans <- x[NA_integer_]
    if (is(x, "Rle"))
        ans <- S4Vectors::Rle(ans)
    ans
}

.zero_cols <- function(x) {
    empty <- new(class(S4Vectors::decode(x)))
    if (is(x, "Rle")) return(S4Vectors::Rle(values = empty))
    empty
}

.na_dframe <- function(dframe, nrows) {
    
    empty <- new("DFrame", nrows = as.integer(nrows))
    
    # IRanges special case for mcols is null
    if (is.null(ncol(dframe))) return(NULL)
    if (ncol(dframe) == 0) return(empty)
    
    tform <- .na_cols
    if (nrows == 0) tform <- .zero_cols
    
    for (i in colnames(dframe)) {
        empty[[i]] <- tform(dframe[[i]])
    }

    empty
    
}

#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectLength
#' @importFrom S4Vectors queryLength
.hits_to_left <- function(x, y, hits, suffix, distance = NULL) {

    left <- x[S4Vectors::subjectHits(hits), ]
    right <- y[S4Vectors::queryHits(hits), ]
    mcols(left) <- .mcols_overlaps_update(left, right, suffix)
    if (!is.null(distance)) mcols(left)$distance <- distance

    # Create empty DF for left entries with non overlap
    only_left <- rep(TRUE, S4Vectors::subjectLength(hits))
    only_left[S4Vectors::subjectHits(hits)] <- FALSE
    rng_only_left <- x[only_left]
    mcols_outer <- .na_dframe(mcols(right), sum(only_left))
    if (!is.null(mcols(rng_only_left))) {
        mcols(rng_only_left) <- cbind(mcols(rng_only_left), mcols_outer)
    } else {
        mcols(rng_only_left) <- mcols_outer
    }
    if (!is.null(distance)) {
        mcols(rng_only_left)$distance <- rep(NA_integer_, sum(only_left))
    }
    names(mcols(rng_only_left)) <- names(mcols(left))

    # Merge left and outer
    left_outer <- c(left, rng_only_left)
    
    # Reorder rows
    left_outer <- left_outer[
        order(c(S4Vectors::subjectHits(hits), which(only_left)))
    ]
    left_outer

}

## Distance between each range of `y` and the anchor of `x` it overlaps, for 
## each pair of `hits` (queries in `y`, subjects in `x`): the given anchor, or 
## the closer of the two. Pairs on different chromosomes (or strands, unless 
## `ignore.strand`) are NA, as in `GenomicRanges::distance()`.
.hits_distance <- function(x, y, hits, anchors, ignore.strand) {
    gr <- y[S4Vectors::queryHits(hits)]
    d <- lapply(anchors, function(i) {
        a <- InteractionSet::anchors(x, type = c("first", "second")[[i]])
        GenomicRanges::distance(
            a[S4Vectors::subjectHits(hits)], gr, ignore.strand = ignore.strand
        )
    })
    do.call(pmin, c(d, na.rm = TRUE))
}
