# Appends distance between interaction anchors

Appends distance between interaction anchors, using
[`InteractionSet::pairdist`](https://rdrr.io/pkg/InteractionSet/man/distances.html)

## Usage

``` r
add_pairdist(x, type = "mid", colname = "pairdist")
```

## Arguments

- x:

  The query GInteractions

- type:

  A character string specifying the type of distance to compute. Can
  take values of "mid", "gap", "span", "diag" or "intra".

- colname:

  name of column to hold pair distance values

## Value

The GInteractions with an additional column containing the distance
between each pair of anchors.

## Examples

``` r
gi <- read.table(text = "
chr1 100 200 chr1 5000 5100 bedpe_example1 30 + -
chr1 1000 5000 chr2 3000 3800 bedpe_example2 100 + -",
col.names = c(
  "seqnames1", "start1", "end1", 
  "seqnames2", "start2", "end2", "name", "score", "strand1", "strand2")
) |> as_ginteractions()

add_pairdist(gi)
#> GInteractions object with 2 interactions and 3 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> |
#>   [1]      chr1   100-200       + ---      chr1 5000-5100       - |
#>   [2]      chr1 1000-5000       + ---      chr2 3000-3800       - |
#>                 name     score  pairdist
#>          <character> <integer> <integer>
#>   [1] bedpe_example1        30      4900
#>   [2] bedpe_example2       100      <NA>
#>   -------
#>   regions: 4 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
