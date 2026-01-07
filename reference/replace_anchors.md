# Replace anchors of a GInteractions

Replace anchors of a GInteractions

## Usage

``` r
replace_anchors(x, id, value)

# S4 method for class 'GInteractions,character,GenomicRanges'
replace_anchors(x, id, value)

# S4 method for class 'GInteractions,numeric,GenomicRanges'
replace_anchors(x, id, value)

# S4 method for class 'PinnedGInteractions,missing,GenomicRanges'
replace_anchors(x, id, value)

# S4 method for class 'AnchoredPinnedGInteractions,missing,GRanges'
replace_anchors(x, id, value)

# S4 method for class 'AnchoredPinnedGInteractions,numeric,GRanges'
replace_anchors(x, id, value)
```

## Arguments

- x:

  a (Pinned)GInteractions object

- id:

  Which anchors to replace ("first" or "second"). Ignored if the
  GInteractions is already pinned to a specific set of anchors.

- value:

  A GRanges object vector the same length as x.

## Value

a (Pinned)GInteractions object.

## Examples

``` r
gi <- read.table(text = "
chr1 11 20 chr1 21 30
chr1 11 20 chr1 51 55
chr1 11 30 chr1 51 55
chr1 11 30 chr2 51 60",
col.names = c(
    "seqnames1", "start1", "end1", 
    "seqnames2", "start2", "end2")
) |> 
  as_ginteractions() |> 
  mutate(type = c('cis', 'cis', 'cis', 'trans'), score = runif(4))

####################################################################
# 1. Replace anchors of a GInteractions object
####################################################################

gi |> replace_anchors(2, value = anchors1(gi))
#> GInteractions object with 4 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     11-20       * ---      chr1     11-20       * |         cis
#>   [2]      chr1     11-20       * ---      chr1     11-20       * |         cis
#>   [3]      chr1     11-30       * ---      chr1     11-30       * |         cis
#>   [4]      chr1     11-30       * ---      chr1     11-30       * |       trans
#>           score
#>       <numeric>
#>   [1]  0.402881
#>   [2]  0.769630
#>   [3]  0.119485
#>   [4]  0.194695
#>   -------
#>   regions: 2 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> replace_anchors(1, value = anchors2(gi))
#> GInteractions object with 4 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     21-30       * ---      chr1     21-30       * |         cis
#>   [2]      chr1     51-55       * ---      chr1     51-55       * |         cis
#>   [3]      chr1     51-55       * ---      chr1     51-55       * |         cis
#>   [4]      chr2     51-60       * ---      chr2     51-60       * |       trans
#>           score
#>       <numeric>
#>   [1]  0.402881
#>   [2]  0.769630
#>   [3]  0.119485
#>   [4]  0.194695
#>   -------
#>   regions: 3 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> replace_anchors(1, value = GenomicRanges::GRanges(c(
  "chr1:1-2", "chr1:2-3", "chr1:3-4", "chr1:4-5"
)))
#> GInteractions object with 4 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1       1-2       * ---      chr1     21-30       * |         cis
#>   [2]      chr1       2-3       * ---      chr1     51-55       * |         cis
#>   [3]      chr1       3-4       * ---      chr1     51-55       * |         cis
#>   [4]      chr1       4-5       * ---      chr2     51-60       * |       trans
#>           score
#>       <numeric>
#>   [1]  0.402881
#>   [2]  0.769630
#>   [3]  0.119485
#>   [4]  0.194695
#>   -------
#>   regions: 7 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

####################################################################
# 2. Replace anchors of a pinned GInteractions object
####################################################################

gi |> pin_by(1) |> replace_anchors(value = anchors1(gi))
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors1
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     11-20       * ---      chr1     21-30       * |         cis
#>   [2]      chr1     11-20       * ---      chr1     51-55       * |         cis
#>   [3]      chr1     11-30       * ---      chr1     51-55       * |         cis
#>   [4]      chr1     11-30       * ---      chr2     51-60       * |       trans
#>           score
#>       <numeric>
#>   [1]  0.402881
#>   [2]  0.769630
#>   [3]  0.119485
#>   [4]  0.194695
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> replace_anchors(1, value = anchors2(gi))
#> GInteractions object with 4 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     21-30       * ---      chr1     21-30       * |         cis
#>   [2]      chr1     51-55       * ---      chr1     51-55       * |         cis
#>   [3]      chr1     51-55       * ---      chr1     51-55       * |         cis
#>   [4]      chr2     51-60       * ---      chr2     51-60       * |       trans
#>           score
#>       <numeric>
#>   [1]  0.402881
#>   [2]  0.769630
#>   [3]  0.119485
#>   [4]  0.194695
#>   -------
#>   regions: 3 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> 
  pin_by(1) |> 
  replace_anchors(value = GenomicRanges::GRanges(c(
    "chr1:1-2", "chr1:2-3", "chr1:3-4", "chr1:4-5"
  ))) |> 
  pin_by(2) |> 
  replace_anchors(value = GenomicRanges::GRanges(c(
    "chr2:1-2", "chr2:2-3", "chr2:3-4", "chr2:4-5"
  ))) 
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors2
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1       1-2       * ---      chr2       1-2       * |         cis
#>   [2]      chr1       2-3       * ---      chr2       2-3       * |         cis
#>   [3]      chr1       3-4       * ---      chr2       3-4       * |         cis
#>   [4]      chr1       4-5       * ---      chr2       4-5       * |       trans
#>           score
#>       <numeric>
#>   [1]  0.402881
#>   [2]  0.769630
#>   [3]  0.119485
#>   [4]  0.194695
#>   -------
#>   regions: 8 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
