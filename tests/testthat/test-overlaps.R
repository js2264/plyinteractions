test_that("plyranges functions work", {
    
    ## Find
    expect_identical(
        find_overlaps(gi, gr) |> length(), 
        5L
    )
    expect_identical(
        find_overlaps_directed(gi, gr) |> length(), 
        3L
    )
    expect_identical(
        gi |> pin_by("first") |> find_overlaps(gr) |> length(), 
        4L
    )
    expect_identical(
        gi |> pin_by("first") |> find_overlaps_directed(gr) |> length(), 
        2L
    )

    ## Count
    expect_identical(
        count_overlaps(gi, gr), 
        c(1L, 1L, 1L, 2L)
    )
    expect_identical(
        count_overlaps_directed(gi, gr), 
        c(1L, 1L, 0L, 1L)
    )
    expect_identical(
        gi |> pin_by("first") |> count_overlaps(gr), 
        c(1L, 1L, 1L, 1L)
    )
    expect_identical(
        gi |> pin_by("first") |> count_overlaps_directed(gr), 
        c(1L, 1L, 0L, 0L)
    )

    ## Filter
    expect_identical(
        filter_by_overlaps(gi, gr) |> length(), 
        4L
    )
    expect_identical(
        filter_by_non_overlaps(gi, gr) |> length(), 
        0L
    )
    expect_identical(
        gi |> pin_by("second") |> filter_by_overlaps(gr) |> length(), 
        2L
    )
    expect_identical(
        gi |> pin_by("second") |> filter_by_non_overlaps(gr) |> length(), 
        2L
    )

    ## Join
    expect_identical(
        join_overlap_left(gi, gr) |> mcols() |> colnames(), 
        c('score', 'type.x', 'id', 'type.y')
    )
    expect_identical(
        join_overlap_left(gi, gr) |> length(), 
        5L
    )
    expect_identical(
        gi |> pin_by(2) |> join_overlap_left(gr) |> mcols() |> colnames(), 
        c('score', 'type.x', 'id', 'type.y')
    )
    expect_identical(
        gi |> pin_by(2) |> join_overlap_left(gr) |> length(), 
        4L
    )
    expect_identical(
        join_overlap_left_directed(gi, gr) |> length(), 
        4L
    )
    expect_identical(
        gi |> pin_by(2) |> join_overlap_left_directed(gr) |> length(), 
        4L
    )

    ## Join, with distance
    expect_identical(
        join_overlap_left(gi, gr, distance = TRUE) |> mcols() |> colnames(), 
        c('score', 'type.x', 'id', 'type.y', 'distance')
    )
    expect_identical(
        join_overlap_left(gi, gr, distance = TRUE)$distance, 
        c(0L, 0L, 0L, 0L, 0L)
    )
    expect_identical(
        join_overlap_left_directed(gi, gr, distance = TRUE)$distance, 
        c(0L, 0L, NA, 0L)
    )
    pinned_join <- gi |> 
        pin_by(2) |> 
        join_overlap_left(gr, maxgap = 25, distance = TRUE)
    expect_identical(
        pinned_join |> mcols() |> colnames(), 
        c('score', 'type.x', 'id', 'type.y', 'distance')
    )
    expect_identical(
        pinned_join$distance, 
        c(0L, 20L, 20L, 0L)
    )
    pinned_join_directed <- gi |> 
        pin_by(2) |> 
        join_overlap_left_directed(gr, maxgap = 25, distance = TRUE)
    expect_identical(
        pinned_join_directed$distance, 
        c(0L, 20L, NA, 0L)
    )

    ## Distance to the closer anchor, or to the pinned one
    gi2 <- InteractionSet::GInteractions(
        GenomicRanges::GRanges("chr1:100-110"), 
        GenomicRanges::GRanges("chr1:150-160")
    )
    gr2 <- GenomicRanges::GRanges("chr1:165-170")
    expect_identical(
        join_overlap_left(gi2, gr2, maxgap = 10, distance = TRUE)$distance, 
        4L
    )
    expect_identical(
        (gi2 |> pin_by(1) |> 
            join_overlap_left(gr2, maxgap = 10, distance = TRUE))$distance, 
        NA_integer_
    )

})
