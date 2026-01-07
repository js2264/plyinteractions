# Find overlaps between a query GInteractions and a GRanges

Find overlaps between a query GInteractions and a GRanges

## Usage

``` r
# S3 method for class 'PinnedGInteractions'
find_overlaps(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y"))

# S3 method for class 'GInteractions'
find_overlaps(x, y, maxgap = -1L, minoverlap = 0L, suffix = c(".x", ".y"))

# S3 method for class 'PinnedGInteractions'
find_overlaps_directed(
  x,
  y,
  maxgap = -1L,
  minoverlap = 0L,
  suffix = c(".x", ".y")
)

# S3 method for class 'GInteractions'
find_overlaps_directed(
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
  `?`[`findOverlaps`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html)
  in the GenomicRanges package for a description of these arguments

- suffix:

  Suffix to add to metadata columns (character vector of length 2,
  default to `c(".x", ".y")`).

## Value

a GInteractions object with rows corresponding to the GInteractions in
`x` that overlap `y`.

## Rationale

`find_overlaps()` will search for any overlap between `GInteractions` in
`x` and `GRanges` in `y`. It will return a `GInteractions` object of
length equal to the number of times `x` overlaps `y`. This
`GInteractions` will have additional metadata columns corresponding to
the metadata from `y`. `find_overlaps_directed()` takes the strandness
of each object into account.

## Pinned `GInteractions`

When using `find_overlaps()` with a `PinnedGInteractions` object, only
the pinned anchors are used to check for overlap with `y`. This is
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
        "seqnames2", "start2", "end2", "strand2"
    )
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
# 1. Find overlaps between GInteractions and a subject GRanges
####################################################################

find_overlaps(gi, gr)
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

find_overlaps_directed(gi, gr)
#> GInteractions object with 1 interaction and 4 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |      id.x
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>   [1]      chr1     11-20       - ---      chr1     21-30       + |         1
#>            type.x      id.y      type.y
#>       <character> <integer> <character>
#>   [1]          gi         1          gr
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

####################################################################
# 2. Find overlaps between PinnedGInteractions and a subject GRanges
####################################################################

gi |> pin_by("first") |> find_overlaps(gr)
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

gi |> pin_by("second") |> find_overlaps(gr)
#> GInteractions object with 2 interactions and 4 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |      id.x
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>   [1]      chr1     11-20       - ---      chr1     21-30       + |         1
#>   [2]      chr1     21-30       - ---      chr2     51-60       + |         4
#>            type.x      id.y      type.y
#>       <character> <integer> <character>
#>   [1]          gi         1          gr
#>   [2]          gi         2          gr
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_by("first") |> find_overlaps_directed(gr)
#> GInteractions object with 0 interactions and 4 metadata columns:
#>    seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |      id.x
#>        <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>         type.x      id.y      type.y
#>    <character> <integer> <character>
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_by("second") |> find_overlaps_directed(gr)
#> GInteractions object with 1 interaction and 4 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |      id.x
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>   [1]      chr1     11-20       - ---      chr1     21-30       + |         1
#>            type.x      id.y      type.y
#>       <character> <integer> <character>
#>   [1]          gi         1          gr
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
