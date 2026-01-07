# Shift pinned anchors of a GInteractions object with plyinteractions

Shift pinned anchors of a GInteractions object with plyinteractions

## Usage

``` r
shift_downstream(x, shift)

# S3 method for class 'Ranges'
shift_downstream(x, shift)

# S3 method for class 'PinnedGInteractions'
shift_downstream(x, shift)

shift_upstream(x, shift)

# S3 method for class 'Ranges'
shift_upstream(x, shift)

# S3 method for class 'PinnedGInteractions'
shift_upstream(x, shift)

shift_right(x, shift)

# S3 method for class 'Ranges'
shift_right(x, shift)

# S3 method for class 'PinnedGInteractions'
shift_right(x, shift)

shift_left(x, shift)

# S3 method for class 'Ranges'
shift_left(x, shift)

# S3 method for class 'PinnedGInteractions'
shift_left(x, shift)
```

## Arguments

- x:

  a PinnedGInteractions object

- shift:

  The amount to move the genomic interval in the Ranges object by.
  Either a non-negative integer vector of length 1 or an integer vector
  the same length as x.

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
# 1. Simple shifting
####################################################################

gi 
#> GInteractions object with 4 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1     11-20       + ---      chr1     21-30       + |  0.996347
#>   [2]      chr1     11-20       + ---      chr1     51-55       + |  0.500192
#>   [3]      chr1     11-30       - ---      chr1     51-55       - |  0.358967
#>   [4]      chr1     11-30       - ---      chr2     51-60       - |  0.774913
#>              type
#>       <character>
#>   [1]         cis
#>   [2]         cis
#>   [3]         cis
#>   [4]       trans
#>   -------
#>   regions: 6 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_by("first") |> shift_left(15)
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors1
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1      -4-5       + ---      chr1     21-30       + |  0.996347
#>   [2]      chr1      -4-5       + ---      chr1     51-55       + |  0.500192
#>   [3]      chr1     -4-15       - ---      chr1     51-55       - |  0.358967
#>   [4]      chr1     -4-15       - ---      chr2     51-60       - |  0.774913
#>              type
#>       <character>
#>   [1]         cis
#>   [2]         cis
#>   [3]         cis
#>   [4]       trans
#>   -------
#>   regions: 6 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_by("second") |> shift_downstream(10)
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors2
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1     11-20       + ---      chr1     31-40       + |  0.996347
#>   [2]      chr1     11-20       + ---      chr1     61-65       + |  0.500192
#>   [3]      chr1     11-30       - ---      chr1     41-45       - |  0.358967
#>   [4]      chr1     11-30       - ---      chr2     41-50       - |  0.774913
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
# 2. Chained shifting of each set of anchors
####################################################################

gi |> 
  pin_by("first") |> shift_downstream(20) |> 
  pin_by("second") |> shift_upstream(20)
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors2
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1     31-40       + ---      chr1      1-10       + |  0.996347
#>   [2]      chr1     31-40       + ---      chr1     31-35       + |  0.500192
#>   [3]      chr1     -9-10       - ---      chr1     71-75       - |  0.358967
#>   [4]      chr1     -9-10       - ---      chr2     71-80       - |  0.774913
#>              type
#>       <character>
#>   [1]         cis
#>   [2]         cis
#>   [3]         cis
#>   [4]       trans
#>   -------
#>   regions: 6 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
