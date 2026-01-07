# Filter GInteractions overlapping with a GRanges

Filter GInteractions overlapping with a GRanges

## Usage

``` r
# S3 method for class 'PinnedGInteractions'
filter_by_overlaps(x, y, maxgap = -1L, minoverlap = 0L)

# S3 method for class 'GInteractions'
filter_by_overlaps(x, y, maxgap = -1L, minoverlap = 0L)

# S3 method for class 'PinnedGInteractions'
filter_by_non_overlaps(x, y, maxgap = -1L, minoverlap = 0L)

# S3 method for class 'GInteractions'
filter_by_non_overlaps(x, y, maxgap = -1L, minoverlap = 0L)
```

## Arguments

- x:

  A (Pinned)GInteractions object

- y:

  A GRanges object

- maxgap, minoverlap:

  See
  `?`[`countOverlaps`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html)
  in the GenomicRanges package for a description of these arguments

## Value

An integer vector of same length as x.

## Pinned `GInteractions`

When using `filter_by_overlaps()` with a `PinnedGInteractions` object,
only the pinned anchors are used to check for overlap with `y`. This is
equivalent to specifying `use.region="both"` in
[`InteractionSet::findOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html).

## Examples

``` r
gi <- read.table(text = "  
    chr1 11 20 - chr1 21 30 + 
    chr1 11 20 - chr1 51 55 + 
    chr1 21 30 - chr1 51 55 + 
    chr1 21 30 - chr2 51 60 +",  
    col.names = c(
        "seqnames1", "start1", "end1", "strand1", 
        "seqnames2", "start2", "end2", "strand2")
) |> as_ginteractions() |> mutate(id = 1:4, type = 'gi')

gr <- GenomicRanges::GRanges(
    c("chr1:20-30:+", "chr2:55-65:-")
) |> mutate(id = 1:2, type = 'gr')

gi
#> GInteractions object with 4 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        id
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>   [1]      chr1     11-20       - ---      chr1     21-30       + |         1
#>   [2]      chr1     11-20       - ---      chr1     51-55       + |         2
#>   [3]      chr1     21-30       - ---      chr1     51-55       + |         3
#>   [4]      chr1     21-30       - ---      chr2     51-60       + |         4
#>              type
#>       <character>
#>   [1]          gi
#>   [2]          gi
#>   [3]          gi
#>   [4]          gi
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gr
#> GRanges object with 2 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        id        type
#>          <Rle> <IRanges>  <Rle> | <integer> <character>
#>   [1]     chr1     20-30      + |         1          gr
#>   [2]     chr2     55-65      - |         2          gr
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

####################################################################
# 1. Filter GInteractions overlapping with a subject GRanges
####################################################################

filter_by_overlaps(gi, gr)
#> GInteractions object with 4 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        id
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>   [1]      chr1     11-20       - ---      chr1     21-30       + |         1
#>   [2]      chr1     11-20       - ---      chr1     51-55       + |         2
#>   [3]      chr1     21-30       - ---      chr1     51-55       + |         3
#>   [4]      chr1     21-30       - ---      chr2     51-60       + |         4
#>              type
#>       <character>
#>   [1]          gi
#>   [2]          gi
#>   [3]          gi
#>   [4]          gi
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

filter_by_non_overlaps(gi, gr)
#> GInteractions object with 0 interactions and 2 metadata columns:
#>    seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        id
#>        <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>           type
#>    <character>
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

####################################################################
# 2. Filter PinnedGInteractions overlapping with a subject GRanges
####################################################################

gi |> pin_by("first") |> filter_by_overlaps(gr)
#> GInteractions object with 4 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        id
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>   [1]      chr1     11-20       - ---      chr1     21-30       + |         1
#>   [2]      chr1     11-20       - ---      chr1     51-55       + |         2
#>   [3]      chr1     21-30       - ---      chr1     51-55       + |         3
#>   [4]      chr1     21-30       - ---      chr2     51-60       + |         4
#>              type
#>       <character>
#>   [1]          gi
#>   [2]          gi
#>   [3]          gi
#>   [4]          gi
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_by("first") |> filter_by_non_overlaps(gr)
#> GInteractions object with 0 interactions and 2 metadata columns:
#>    seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        id
#>        <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>           type
#>    <character>
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_by("second") |> filter_by_overlaps(gr)
#> GInteractions object with 2 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        id
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>   [1]      chr1     11-20       - ---      chr1     21-30       + |         1
#>   [2]      chr1     21-30       - ---      chr2     51-60       + |         4
#>              type
#>       <character>
#>   [1]          gi
#>   [2]          gi
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_by("second") |> filter_by_non_overlaps(gr)
#> GInteractions object with 2 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        id
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>   [1]      chr1     11-20       - ---      chr1     51-55       + |         2
#>   [2]      chr1     21-30       - ---      chr1     51-55       + |         3
#>              type
#>       <character>
#>   [1]          gi
#>   [2]          gi
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
