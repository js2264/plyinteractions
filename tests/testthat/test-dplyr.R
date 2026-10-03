test_that("dplyr functions work", {
    
    ## arrange

    expect_identical(
        gi |> mutate(idx = seq(1, length(gi))) |> 
            arrange(score) |> 
            strand1(), 
        S4Vectors::Rle(factor(c('+', '-', '+', '-'), c('+', '-', '*')))
    )

    ## tally/count

    ggi |> count(strand1) |> expect_no_error()
    expect_identical(
        ggi |> count(group), 
        new("DFrame", rownames = NULL, nrows = 2L, elementType = "ANY", 
            elementMetadata = NULL, metadata = list(), listData = list(
                group = c(1, 2), n = c(2L, 2L)))
    )
    expect_identical(
        ggi |> count(type), 
        new("DFrame", rownames = NULL, nrows = 3L, elementType = "ANY", 
            elementMetadata = NULL, metadata = list(), listData = list(
                group = c(1, 2, 2), type = c("cis", "cis", "trans"), 
                n = c(2L, 1L, 1L)))
    )
    expect_identical(
        ggi |> tally(), 
        new("DFrame", rownames = NULL, nrows = 2L, elementType = "ANY", 
            elementMetadata = NULL, metadata = list(), listData = list(
                group = c(1, 2), n = c(2L, 2L)))
    )
    expect_identical(
        ggi |> count(strand1), 
        new("DFrame", rownames = NULL, nrows = 2L, elementType = "ANY", 
            elementMetadata = NULL, metadata = list(), listData = list(
                group = c(1, 2), strand1 = new("Rle", values = structure(1:2, levels = c("+", 
                "-", "*"), class = "factor"), lengths = c(1L, 1L), elementMetadata = NULL, 
                    metadata = list()), n = c(2L, 2L)))
    )
    expect_identical(
        gi |> count(strand1), 
        new("DFrame", rownames = NULL, nrows = 2L, elementType = "ANY", 
            elementMetadata = NULL, metadata = list(), listData = list(
                strand1 = new("Rle", values = structure(1:2, levels = c("+", 
                "-", "*"), class = "factor"), lengths = c(1L, 1L), elementMetadata = NULL, 
                    metadata = list()), n = c(2L, 2L)))
    )
    expect_identical(
        ggi |> count(strand1, wt = score), 
        new("DFrame", rownames = NULL, nrows = 2L, elementType = "ANY", 
            elementMetadata = NULL, metadata = list(), listData = list(
                group = c(1, 2), strand1 = new("Rle", values = structure(1:2, levels = c("+", 
                "-", "*"), class = "factor"), lengths = c(1L, 1L), elementMetadata = NULL, 
                    metadata = list()), n = c(0.736002816120163, 1.23265417455696
                ))), 
        tolerance = 1e-4
    )
    expect_identical(
        gi |> count(strand1, wt = score), 
        new("DFrame", rownames = NULL, nrows = 2L, elementType = "ANY", 
            elementMetadata = NULL, metadata = list(), listData = list(
                strand1 = new("Rle", values = structure(1:2, levels = c("+", 
                "-", "*"), class = "factor"), lengths = c(1L, 1L), elementMetadata = NULL, 
                    metadata = list()), n = c(0.736002816120163, 1.23265417455696
                ))), 
        tolerance = 1e-4
    )
    expect_identical(
        ggi |> count(strand1, wt = score, sort = TRUE), 
        new("DFrame", rownames = NULL, nrows = 2L, elementType = "ANY", 
            elementMetadata = NULL, metadata = list(), listData = list(
                group = c(2, 1), strand1 = new("Rle", values = structure(2:1, levels = c("+", 
                "-", "*"), class = "factor"), lengths = c(1L, 1L), elementMetadata = NULL, 
                    metadata = list()), n = c(1.23265417455696, 0.736002816120163
                ))), 
        tolerance = 1e-4
    )
    expect_identical(
        gi |> count(strand1, wt = score, sort = TRUE), 
        new("DFrame", rownames = NULL, nrows = 2L, elementType = "ANY", 
            elementMetadata = NULL, metadata = list(), listData = list(
                strand1 = new("Rle", values = structure(2:1, levels = c("+", 
                "-", "*"), class = "factor"), lengths = c(1L, 1L), elementMetadata = NULL, 
                    metadata = list()), n = c(1.23265417455696, 0.736002816120163
                ))), 
        tolerance = 1e-4
    )

    ## filter

    expect_identical(
        gi |> filter(strand1 == '-') |> ranges1(), 
        new("IRanges", start = c(11L, 11L), width = c(20L, 20L), NAMES = NULL, 
            elementType = "ANY", elementMetadata = NULL, metadata = list())
    )
    expect_identical(
        gi |> filter(strand1 == '+') |> ranges1(), 
        new("IRanges", start = c(11L, 11L), width = c(10L, 10L), NAMES = NULL, 
            elementType = "ANY", elementMetadata = NULL, metadata = list())
    )
    ggi |> filter(strand1 == '+') |> expect_error()

    ## mutate

    expect_identical(
        gi |> mutate(strand1 = '-') |> strand1(), 
        new("Rle", values = structure(2L, levels = c("+", "-", "*"), class = "factor"), 
        lengths = 4L, elementMetadata = NULL, metadata = list())
    )
    expect_identical(
        gi |> mutate(xxx = 1) |> mcols() |> subset(, c(2, 3)), 
        new("DFrame", rownames = NULL, nrows = 4L, elementType = "ANY", 
        elementMetadata = NULL, metadata = list(), listData = list(
            type = c("cis", "cis", "cis", "trans"), xxx = c(1, 1, 
            1, 1)))
    )
    expect_identical(
        mutate(gi, xxx = IRanges::RleList(c(1, 2), c(3, 4)))$xxx,
        IRanges::RleList(c(1, 2), c(3, 4), c(1, 2), c(3, 4))
    )

    ## mutate, on grouped GInteractions: evaluated within each group, and 
    ## grouped again
    ggi_mutated <- ggi |> mutate(m = mean(score), m2 = m * 2)
    expect_s4_class(ggi_mutated, "GroupedGInteractions")
    expect_identical(group_vars(ggi_mutated), "group")
    expect_equal(
        ggi_mutated$m, 
        rep(c(mean(gi$score[1:2]), mean(gi$score[3:4])), each = 2)
    )
    expect_equal(ggi_mutated$m2, ggi_mutated$m * 2)
    expect_identical(
        ggi |> mutate(strand1 = '-') |> strand1() |> as.character(), 
        rep('-', 4)
    )
    expect_s4_class(ggi |> mutate(strand1 = '-'), "GroupedGInteractions")
    expect_identical(ggi |> mutate(group = 1) |> n_groups(), 1L)

    ## mutate, on pinned GInteractions: they stay pinned (and anchored)
    pgi_mutated <- pgi |> mutate(s2 = score * 2, start2 = 1)
    expect_s4_class(pgi_mutated, "PinnedGInteractions")
    expect_identical(pin(pgi_mutated), 2L)
    expect_identical(pgi_mutated$s2, gi$score * 2)
    expect_identical(start2(pgi_mutated), rep(1L, 4))
    apgi_mutated <- apgi |> mutate(s2 = score * 2, start1 = 1)
    expect_s4_class(apgi_mutated, "AnchoredPinnedGInteractions")
    expect_identical(anchor(apgi_mutated), "5p")
    expect_identical(apgi_mutated$s2, gi$score * 2)
    expect_identical(start1(apgi_mutated), rep(1L, 4))

    ## subsetting, on pinned GInteractions: they stay pinned (and anchored)
    expect_s4_class(pgi[2:3], "PinnedGInteractions")
    expect_identical(unpin(pgi[2:3]), gi[2:3])
    expect_identical(pin(filter(pgi, score > 0.5)), 2L)
    expect_identical(
        pgi |> filter(score > 0.5) |> unpin(), 
        gi |> filter(score > 0.5)
    )
    expect_identical(apgi |> slice(2:3) |> anchor(), "5p")
    expect_identical(apgi |> slice(2:3) |> unpin(), gi |> slice(2:3))
    expect_identical(pgi |> arrange(score) |> unpin(), gi |> arrange(score))
    expect_identical(
        apgi |> mutate(width2 = 100) |> width2(), 
        c(100L, 100L, 100L, 100L)
    )
    expect_s4_class(
        apgi |> mutate(width2 = 100), 
        "AnchoredPinnedGInteractions"
    )

    ## rename

    expect_error(gi |> rename(strand1 = strand2))
    expect_error(gi |> rename(yyy = xxx))
    expect_identical(
        gi |> rename(xx = type) |> as_tibble() |> colnames(), 
        c("seqnames1", "start1", "end1", "width1", "strand1", "seqnames2", 
        "start2", "end2", "width2", "strand2", "score", "xx")
    )

    ## as_tibble, length and metadata columns of grouped and pinned 
    ## GInteractions
    expect_identical(
        as_tibble(ggi), 
        as_tibble(ungroup(ggi))
    )
    expect_identical(length(ggi), 4L)
    expect_identical(length(pgi), 4L)
    expect_identical(length(apgi), 4L)
    pgi2 <- pgi
    pgi2$new <- 1
    expect_s4_class(pgi2, "PinnedGInteractions")
    expect_identical(pgi2$new, rep(1, 4))
    ggi2 <- ggi
    ggi2$group <- c(1, 2, 3, 4)
    expect_s4_class(ggi2, "GroupedGInteractions")
    expect_identical(n_groups(ggi2), 4L)
    expect_identical(
        as_tibble(pgi), 
        as_tibble(gi)
    )
    expect_identical(
        as_tibble(apgi), 
        as_tibble(gi)
    )
    expect_identical(
        dim(as_tibble(ggi)), 
        c(4L, 13L)
    )

    ## select
    expect_identical(
        gi |> select(type) |> as_tibble() |> colnames(), 
        c("seqnames1", "start1", "end1", "width1", "strand1", "seqnames2", 
        "start2", "end2", "width2", "strand2", "type")
    )
    expect_identical(
        gi |> select(type, .drop_ranges = TRUE) |> as_tibble() |> colnames(), 
        c("type")
    )
    expect_error(gi |> select(strand1))
    expect_identical(
        gi |> select(strand1, .drop_ranges = TRUE) |> as_tibble() |> colnames(), 
        c("strand1")
    )

    ## slice
    expect_identical(
        gi |> slice(1:3) |> start1(), 
        c(11L, 11L, 11L)
    )
    gi |> slice(1:5) |> expect_error()
    gi |> slice('error') |> expect_error()

    ## summarize
    expect_identical(
        ggi |> summarize(m = mean(score)), 
        new("DFrame", rownames = NULL, nrows = 2L, elementType = "ANY", 
        elementMetadata = NULL, metadata = list(), listData = list(
            group = c(1, 2), m = c(0.368001408060081, 0.616327087278478
            ))), 
        tolerance = 1e-4
    )
    expect_identical(
        ggi |> summarize(m = table(type)), 
        new("DFrame", rownames = NULL, nrows = 2L, elementType = "ANY", 
            elementMetadata = NULL, metadata = list(), listData = list(
                group = c(1, 2), m = new("CompressedIntegerList", elementType = "integer", 
                    elementMetadata = NULL, metadata = list(), unlistData = c(cis = 2L, 
                    trans = 0L, cis = 1L, trans = 1L), partitioning = new("PartitioningByEnd", 
                        end = c(2L, 4L), NAMES = c("1", "2"), elementType = "ANY", 
                        elementMetadata = NULL, metadata = list()))))
    )
    expect_identical(
        ggi |> summarize(m = table(type), n = table(group)), 
        new("DFrame", rownames = NULL, nrows = 2L, elementType = "ANY", 
            elementMetadata = NULL, metadata = list(), listData = list(
                group = c(1, 2), m = new("CompressedIntegerList", elementType = "integer", 
                    elementMetadata = NULL, metadata = list(), unlistData = c(cis = 2L, 
                    trans = 0L, cis = 1L, trans = 1L), partitioning = new("PartitioningByEnd", 
                        end = c(2L, 4L), NAMES = c("1", "2"), elementType = "ANY", 
                        elementMetadata = NULL, metadata = list())), 
                n = new("CompressedIntegerList", elementType = "integer", 
                    elementMetadata = NULL, metadata = list(), unlistData = c(`1` = 2L, 
                    `2` = 0L, `1` = 0L, `2` = 2L), partitioning = new("PartitioningByEnd", 
                        end = c(2L, 4L), NAMES = c("1", "2"), elementType = "ANY", 
                        elementMetadata = NULL, metadata = list()))))
    )

})

test_that("dplyr group functions work", {

    expect_identical(
        group_data(ggi), 
        new("DFrame", rownames = NULL, nrows = 2L, elementType = "ANY", 
        elementMetadata = NULL, metadata = list(), listData = list(
            group = c(1, 2), .rows = new("CompressedIntegerList", 
                elementType = "integer", elementMetadata = NULL, 
                metadata = list(), unlistData = 1:4, partitioning = new("PartitioningByEnd", 
                    end = c(2L, 4L), NAMES = NULL, elementType = "ANY", 
                    elementMetadata = NULL, metadata = list()))))
    )
    expect_identical(
        group_keys(ggi), 
        new("DFrame", rownames = NULL, nrows = 2L, elementType = "ANY", 
            elementMetadata = NULL, metadata = list(), listData = list(
                group = c(1, 2)))
    )
    expect_identical(
        group_indices(ggi), 
        S4Vectors::Rle(c(1L, 1L, 2L, 2L))
    )
    expect_identical(
        group_vars(ggi), 
        "group"
    )
    expect_identical(
        groups(ggi), 
        list(rlang::sym("group"))
    )
    expect_identical(
        group_size(ggi), 
        c(2L, 2L)
    )
    expect_identical(
        n_groups(ggi), 
        2L
    )

})
