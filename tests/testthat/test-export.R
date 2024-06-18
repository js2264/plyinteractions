test_that("export functions work", {
    
    ## write_bedpe
    gi |> write_bedpe(file = 'gi.bedpe') |> expect_true()
    gi |> write_bedpe(file = 'gi.bedpe', scores = 'score') |> expect_true()
    gi |> write_bedpe(file = 'gi.bedpe', scores = 'fail') |> expect_error()

    ## write_pairs
    gi |> write_pairs(file = 'gi.pairs') |> expect_message()
    gi |>
        write_pairs(file = 'gi.pairs', seqlengths = c(I = 100, II = 100)) |> 
        expect_error()
    gi |>
        write_pairs(file = 'gi.pairs', seqlengths = c(chr1 = 100, chr2 = 100)) |> 
        expect_true()
    gi |>
        write_pairs(file = 'gi.pairs', seqlengths = c(chr1 = 100, chr3 = 30)) |> 
        expect_error()
    gi |>
        write_pairs(file = 'gi.pairs', seqlengths = c(chr1 = 100, chr2 = 100, chr3 = 30)) |> 
        expect_true()
    gi |>
        write_pairs(file = 'gi.pairs', seqlengths = c(chr1 = 100, chr2 = 30)) |> 
        expect_error()

})
