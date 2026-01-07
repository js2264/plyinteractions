# Join overlaps between a query GInteractions and a GRanges

Join overlaps between a query GInteractions and a GRanges

## Usage

``` r
# S3 method for class 'PinnedGInteractions'
join_overlap_left(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y"))

# S3 method for class 'GInteractions'
join_overlap_left(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y"))

# S3 method for class 'PinnedGInteractions'
join_overlap_left_directed(
  x,
  y,
  maxgap = -1L,
  minoverlap = 0L,
  suffix = c(".x", ".y")
)

# S3 method for class 'GInteractions'
join_overlap_left_directed(
  x,
  y,
  maxgap = -1L,
  minoverlap = 0L,
  suffix = c(".x", ".y")
)
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

- suffix:

  Suffix to add to metadata columns (character vector of length 2,
  default to `c(".x", ".y")`).

## Value

An integer vector of same length as x.

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
# 1. Join overlaps between GInteractions and a subject GRanges
####################################################################

join_overlap_left(gi, gr)
#> GInteractions object with 5 interactions and 4 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |      id.x
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>   [1]      chr1     11-20       - ---      chr1     21-30       + |         1
#>   [2]      chr1     11-20       - ---      chr1     51-55       + |         2
#>   [3]      chr1     21-30       - ---      chr1     51-55       + |         3
#>   [4]      chr1     21-30       - ---      chr2     51-60       + |         4
#>   [5]      chr1     21-30       - ---      chr2     51-60       + |         4
#>            type.x      id.y      type.y
#>       <character> <integer> <character>
#>   [1]          gi         1          gr
#>   [2]          gi         1          gr
#>   [3]          gi         1          gr
#>   [4]          gi         1          gr
#>   [5]          gi         2          gr
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

join_overlap_left_directed(gi, gr)
#> GInteractions object with 4 interactions and 4 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |      id.x
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>   [1]      chr1     11-20       - ---      chr1     21-30       + |         1
#>   [2]      chr1     11-20       - ---      chr1     51-55       + |         2
#>   [3]      chr1     21-30       - ---      chr1     51-55       + |         3
#>   [4]      chr1     21-30       - ---      chr2     51-60       + |         4
#>            type.x      id.y      type.y
#>       <character> <integer> <character>
#>   [1]          gi         1          gr
#>   [2]          gi      <NA>        <NA>
#>   [3]          gi      <NA>        <NA>
#>   [4]          gi      <NA>        <NA>
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

####################################################################
# 2. Join overlaps between PinnedGInteractions and a subject GRanges
####################################################################

gi |> pin_by("first") |> join_overlap_left(gr)
#> GInteractions object with 4 interactions and 4 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |      id.x
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>   [1]      chr1     11-20       - ---      chr1     21-30       + |         1
#>   [2]      chr1     11-20       - ---      chr1     51-55       + |         2
#>   [3]      chr1     21-30       - ---      chr1     51-55       + |         3
#>   [4]      chr1     21-30       - ---      chr2     51-60       + |         4
#>            type.x      id.y      type.y
#>       <character> <integer> <character>
#>   [1]          gi         1          gr
#>   [2]          gi         1          gr
#>   [3]          gi         1          gr
#>   [4]          gi         1          gr
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_by("first") |> join_overlap_left_directed(gr)
#> GInteractions object with 4 interactions and 4 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |      id.x
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>   [1]      chr1     11-20       - ---      chr1     21-30       + |         1
#>   [2]      chr1     11-20       - ---      chr1     51-55       + |         2
#>   [3]      chr1     21-30       - ---      chr1     51-55       + |         3
#>   [4]      chr1     21-30       - ---      chr2     51-60       + |         4
#>            type.x      id.y      type.y
#>       <character> <integer> <character>
#>   [1]          gi      <NA>        <NA>
#>   [2]          gi      <NA>        <NA>
#>   [3]          gi      <NA>        <NA>
#>   [4]          gi      <NA>        <NA>
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_by("second") |> join_overlap_left(gr)
#> GInteractions object with 4 interactions and 4 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |      id.x
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>   [1]      chr1     11-20       - ---      chr1     21-30       + |         1
#>   [2]      chr1     11-20       - ---      chr1     51-55       + |         2
#>   [3]      chr1     21-30       - ---      chr1     51-55       + |         3
#>   [4]      chr1     21-30       - ---      chr2     51-60       + |         4
#>            type.x      id.y      type.y
#>       <character> <integer> <character>
#>   [1]          gi         1          gr
#>   [2]          gi      <NA>        <NA>
#>   [3]          gi      <NA>        <NA>
#>   [4]          gi         2          gr
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_by("second") |> join_overlap_left_directed(gr)
#> GInteractions object with 4 interactions and 4 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |      id.x
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>   [1]      chr1     11-20       - ---      chr1     21-30       + |         1
#>   [2]      chr1     11-20       - ---      chr1     51-55       + |         2
#>   [3]      chr1     21-30       - ---      chr1     51-55       + |         3
#>   [4]      chr1     21-30       - ---      chr2     51-60       + |         4
#>            type.x      id.y      type.y
#>       <character> <integer> <character>
#>   [1]          gi         1          gr
#>   [2]          gi      <NA>        <NA>
#>   [3]          gi      <NA>        <NA>
#>   [4]          gi      <NA>        <NA>
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
