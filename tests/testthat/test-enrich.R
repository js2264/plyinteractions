test_that("enrich functions work", {
    
    gi <- read.table(text = "
        chr1 100 200 chr1 5000 5100 bedpe_example1 30 + -
        chr1 1000 5000 chr2 3000 3800 bedpe_example2 100 + -",
        col.names = c(
            "seqnames1", "start1", "end1", 
            "seqnames2", "start2", "end2", "name", "score", "strand1", "strand2"
        )
    ) |> as_ginteractions()

    add_pairdist(gi) |> expect_identical(
        new("GInteractions", anchor1 = 1:2, anchor2 = 3:4, regions = new("GRanges", 
            seqnames = new("Rle", values = structure(1:2, levels = c("chr1", 
            "chr2"), class = "factor"), lengths = c(3L, 1L), elementMetadata = NULL, 
                metadata = list()), ranges = new("IRanges", start = c(100L, 
            1000L, 5000L, 3000L), width = c(101L, 4001L, 101L, 801L), 
                NAMES = NULL, elementType = "ANY", elementMetadata = NULL, 
                metadata = list()), strand = new("Rle", values = structure(1:2, levels = c("+", 
            "-", "*"), class = "factor"), lengths = c(2L, 2L), elementMetadata = NULL, 
                metadata = list()), seqinfo = new("Seqinfo", seqnames = c("chr1", 
            "chr2"), seqlengths = c(NA_integer_, NA_integer_), is_circular = c(NA, 
            NA), genome = c(NA_character_, NA_character_)), elementMetadata = new("DFrame", 
                rownames = NULL, nrows = 4L, elementType = "ANY", elementMetadata = NULL, 
                metadata = list(), listData = structure(list(), names = character(0))), 
            elementType = "ANY", metadata = list()), NAMES = NULL, elementMetadata = new("DFrame", 
            rownames = NULL, nrows = 2L, elementType = "ANY", elementMetadata = NULL, 
            metadata = list(), listData = list(name = c("bedpe_example1", 
            "bedpe_example2"), score = c(30L, 100L), pairdist = c(4900L, 
            NA))), metadata = list())
        )
    
    add_pairdist(gi, colname = 's') |> expect_identical(
        new("GInteractions", anchor1 = 1:2, anchor2 = 3:4, regions = new("GRanges", 
            seqnames = new("Rle", values = structure(1:2, levels = c("chr1", 
            "chr2"), class = "factor"), lengths = c(3L, 1L), elementMetadata = NULL, 
                metadata = list()), ranges = new("IRanges", start = c(100L, 
            1000L, 5000L, 3000L), width = c(101L, 4001L, 101L, 801L), 
                NAMES = NULL, elementType = "ANY", elementMetadata = NULL, 
                metadata = list()), strand = new("Rle", values = structure(1:2, levels = c("+", 
            "-", "*"), class = "factor"), lengths = c(2L, 2L), elementMetadata = NULL, 
                metadata = list()), seqinfo = new("Seqinfo", seqnames = c("chr1", 
            "chr2"), seqlengths = c(NA_integer_, NA_integer_), is_circular = c(NA, 
            NA), genome = c(NA_character_, NA_character_)), elementMetadata = new("DFrame", 
                rownames = NULL, nrows = 4L, elementType = "ANY", elementMetadata = NULL, 
                metadata = list(), listData = structure(list(), names = character(0))), 
            elementType = "ANY", metadata = list()), NAMES = NULL, elementMetadata = new("DFrame", 
            rownames = NULL, nrows = 2L, elementType = "ANY", elementMetadata = NULL, 
            metadata = list(), listData = list(name = c("bedpe_example1", 
            "bedpe_example2"), score = c(30L, 100L), s = c(4900L, 
            NA))), metadata = list())
        )
    
    gr <- read.table(text = "
        chr1 100 200
        chr1 5000 5100
        chr1 1000 5000
        chr2 3000 3800",
        col.names = c(
        "seqnames", "start", "end" 
    )) |> plyranges::as_granges()

    pair_granges(gr) |> expect_identical(
        new("GInteractions", anchor1 = c(1L, 1L, 1L, 3L, 3L, 2L), anchor2 = c(3L, 
            2L, 4L, 2L, 4L, 4L), regions = new("GRanges", seqnames = new("Rle", 
            values = structure(1:2, levels = c("chr1", "chr2"), class = "factor"), 
            lengths = c(3L, 1L), elementMetadata = NULL, metadata = list()), 
            ranges = new("IRanges", start = c(100L, 1000L, 5000L, 3000L
            ), width = c(101L, 4001L, 101L, 801L), NAMES = NULL, elementType = "ANY", 
                elementMetadata = NULL, metadata = list()), strand = new("Rle", 
                values = structure(3L, levels = c("+", "-", "*"), class = "factor"), 
                lengths = 4L, elementMetadata = NULL, metadata = list()), 
            seqinfo = new("Seqinfo", seqnames = c("chr1", "chr2"), seqlengths = c(NA_integer_, 
            NA_integer_), is_circular = c(NA, NA), genome = c(NA_character_, 
            NA_character_)), elementMetadata = new("DFrame", rownames = NULL, 
                nrows = 4L, elementType = "ANY", elementMetadata = NULL, 
                metadata = list(), listData = structure(list(), names = character(0))), 
            elementType = "ANY", metadata = list()), NAMES = NULL, elementMetadata = new("DFrame", 
            rownames = NULL, nrows = 6L, elementType = "ANY", elementMetadata = NULL, 
            metadata = list(), listData = structure(list(), names = character(0))), 
            metadata = list()
        )
    )
})
