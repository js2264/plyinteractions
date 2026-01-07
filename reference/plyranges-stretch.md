# Stretch pinned anchors of a GInteractions object with plyranges

Stretch pinned anchors of a GInteractions object with plyranges

## Usage

``` r
# S3 method for class 'AnchoredPinnedGInteractions'
stretch(x, extend)

# S3 method for class 'PinnedGInteractions'
stretch(x, extend)
```

## Arguments

- x:

  a PinnedGInteractions object

- extend:

  The amount to alter the width of a Ranges object by. Either an integer
  vector of length 1 or an integer vector the same length as x.

## Value

A PinnedGInteractions object

## Examples

``` r
gi <- read.table(text = "
chr1 11 20 chr1 21 30 + +
chr1 11 20 chr1 51 55 + +
chr1 11 30 chr1 51 55 - -
chr1 11 30 chr2 51 60 - -",
col.names = c(
  "seqnames1", "start1", "end1", 
  "seqnames2", "start2", "end2", "strand1", "strand2")
) |> 
  as_ginteractions() |> 
  mutate(score = runif(4), type = c('cis', 'cis', 'cis', 'trans'))

####################################################################
# 1. Simple stretching
####################################################################

gi 
#> GInteractions object with 4 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1     11-20       + ---      chr1     21-30       + |  0.858666
#>   [2]      chr1     11-20       + ---      chr1     51-55       + |  0.566894
#>   [3]      chr1     11-30       - ---      chr1     51-55       - |  0.252997
#>   [4]      chr1     11-30       - ---      chr2     51-60       - |  0.918803
#>              type
#>       <character>
#>   [1]         cis
#>   [2]         cis
#>   [3]         cis
#>   [4]       trans
#>   -------
#>   regions: 6 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_by("first") |> anchor_start() |> stretch(15)
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors1
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1     11-35       + ---      chr1     21-30       + |  0.858666
#>   [2]      chr1     11-35       + ---      chr1     51-55       + |  0.566894
#>   [3]      chr1     11-45       - ---      chr1     51-55       - |  0.252997
#>   [4]      chr1     11-45       - ---      chr2     51-60       - |  0.918803
#>              type
#>       <character>
#>   [1]         cis
#>   [2]         cis
#>   [3]         cis
#>   [4]       trans
#>   -------
#>   regions: 6 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_by("second") |> anchor_center() |> stretch(10)
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors2
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1     11-20       + ---      chr1     16-35       + |  0.858666
#>   [2]      chr1     11-20       + ---      chr1     46-60       + |  0.566894
#>   [3]      chr1     11-30       - ---      chr1     46-60       - |  0.252997
#>   [4]      chr1     11-30       - ---      chr2     46-65       - |  0.918803
#>              type
#>       <character>
#>   [1]         cis
#>   [2]         cis
#>   [3]         cis
#>   [4]       trans
#>   -------
#>   regions: 6 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_by("second") |> anchor_3p() |> stretch(20)
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors2
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1     11-20       + ---      chr1      1-30       + |  0.858666
#>   [2]      chr1     11-20       + ---      chr1     31-55       + |  0.566894
#>   [3]      chr1     11-30       - ---      chr1     51-75       - |  0.252997
#>   [4]      chr1     11-30       - ---      chr2     51-80       - |  0.918803
#>              type
#>       <character>
#>   [1]         cis
#>   [2]         cis
#>   [3]         cis
#>   [4]       trans
#>   -------
#>   regions: 6 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

####################################################################
# 2. Chained stretching of each set of anchors
####################################################################

gi |> 
  pin_by("first") |> anchor_start() |> stretch(20) |> 
  pin_by("second") |> stretch(20)
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors2
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1     11-40       + ---      chr1     11-40       + |  0.858666
#>   [2]      chr1     11-40       + ---      chr1     41-65       + |  0.566894
#>   [3]      chr1     11-50       - ---      chr1     41-65       - |  0.252997
#>   [4]      chr1     11-50       - ---      chr2     41-70       - |  0.918803
#>              type
#>       <character>
#>   [1]         cis
#>   [2]         cis
#>   [3]         cis
#>   [4]       trans
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
