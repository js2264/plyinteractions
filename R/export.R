#' Export GInteractions as `bedpe` or `pairs` files
#'
#' @description `write_*` functions are provided to export a GInteractions 
#' object into these two file formats. See 4DN documentation 
#' (https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md)
#' and UCSC documentation 
#' (https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format)
#' for more details. 
#'
#' @param x a GInteractions object.
#' @param file path to a `.bedpe` or `.pairs` file to save the genomic interactions.
#'
#' @return TRUE
#'
#' @rdname ginteractions-export
#' 
#' @examples
#' gi <- read.table(text = "
#' chr1 100 200 chr1 5000 5100 bedpe_example1 30 + -
#' chr1 1000 5000 chr1 3000 3800 bedpe_example2 100 + -",
#' col.names = c(
#'   "seqnames1", "start1", "end1", 
#'   "seqnames2", "start2", "end2", "name", "score", "strand1", "strand2")
#' ) |> as_ginteractions()
#' 
#' write_bedpe(gi, 'gi.bedpe')
#' write_pairs(gi, 'gi.pairs')

write_bedpe <- function(x, file, scores = NULL) {
        
    tab <- x |> 
        tibble::as_tibble() |> 
        dplyr::select(seqnames1, start1, end1, seqnames2, start2, end2, strand1, strand2) |> 
        dplyr::mutate(name = ".", score = ".") |> 
        dplyr::relocate(name, score, .after = end2) |> 
        dplyr::mutate(start1 = start1-1, start2 = start2-1)
    
    if (!is.null(scores)) {
        if (scores %in% colnames(GenomicRanges::mcols(x))) {
            tab[[scores]] <- GenomicRanges::mcols(x)[[scores]]
        } 
        else {
            stop("column `", scores, "` does not exist.")
        }
    }

    cn <- colnames(GenomicRanges::mcols(x))
    cn <- cn[!cn %in% colnames(tab)]
    d <- GenomicRanges::mcols(x)[,cn] |> as.data.frame()
    colnames(d) <- cn
    tab <- cbind(tab, d)

    utils::write.table(
        tab, 
        file = file, quote = FALSE, sep = '\t', 
        row.names = FALSE, col.names = FALSE
    )

    TRUE
}

write_pairs <- function(x, file, seqlengths = GenomeInfoDb::seqlengths(x), scores = NULL) {
    
    if (any(is.na(seqlengths))) {
        message("No `seqlengths` provided. The `chromsizes` will be inferred from the interactions and will most likely be inaccurate.")
        seqlengths <- InteractionSet::regions(x) |> 
            as_tibble() |> 
            dplyr::group_by(seqnames) |> 
            dplyr::select(seqnames, end) |> 
            dplyr::slice_max(end, n = 1, with_ties = FALSE) |> 
            tibble::deframe()
    }
    else {
        seqlengths_detected <- InteractionSet::regions(x) |> 
            as_tibble() |> 
            dplyr::group_by(seqnames) |> 
            dplyr::select(seqnames, end) |> 
            dplyr::slice_max(end, n = 1, with_ties = FALSE) |> 
            tibble::deframe()
        if (any(!names(seqlengths_detected) %in% names(seqlengths))) {
            stop("Missing seqlengths in the manually provided `seqlengths`")
        }
        if (any(!(seqlengths[names(seqlengths_detected)] >= seqlengths_detected))) {
            cat("Detected `seqlengths:`\n")
            print(seqlengths_detected)
            cat("Provided `seqlengths:`\n")
            print(seqlengths)
            stop("Interactions exceed provided seqlengths")
        }
    }
    chromsizes <- lapply(seq_along(seqlengths), function(.i) {
        .x <- seqlengths[.i]
        .y <- names(seqlengths)[.i]
        stringr::str_glue("#chromsize: {.y} {.x}")
    }) |> unlist()

    header <- c("## pairs format v1.0", chromsizes, "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2")
    writeLines(header, file)

    tab <- x |> 
        tibble::as_tibble() |> 
        dplyr::select(seqnames1, start1, seqnames2, start2, strand1, strand2) |> 
        dplyr::mutate(start1 = start1-1, start2 = start2-1) |> 
        dplyr::mutate(name = dplyr::n() |> seq_len()) |> 
        dplyr::relocate(name, before = seqnames1)
    
    utils::write.table(
        tab, 
        file = file, quote = FALSE, sep = ' ', 
        row.names = FALSE, col.names = FALSE, 
        append = TRUE
    )

    TRUE
}
