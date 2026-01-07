# Generate flanking regions from pinned anchors of a GInteractions object

Generate flanking regions from pinned anchors of a GInteractions object

## Usage

``` r
flank_downstream(x, width)

# S3 method for class 'Ranges'
flank_downstream(x, width)

# S3 method for class 'PinnedGInteractions'
flank_downstream(x, width)

flank_left(x, width)

# S3 method for class 'Ranges'
flank_left(x, width)

# S3 method for class 'PinnedGInteractions'
flank_left(x, width)

flank_upstream(x, width)

# S3 method for class 'Ranges'
flank_upstream(x, width)

# S3 method for class 'PinnedGInteractions'
flank_upstream(x, width)

flank_right(x, width)

# S3 method for class 'Ranges'
flank_right(x, width)

# S3 method for class 'PinnedGInteractions'
flank_right(x, width)
```

## Arguments

- x:

  a PinnedGInteractions object

- width:

  The width of the flanking region relative to the ranges in x. Either
  an integer vector of length 1 or an integer vector the same length
  as x. The width can be negative in which case the flanking region is
  reversed.

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
# 1. Simple flanking
####################################################################

gi 
#> GInteractions object with 4 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1     11-20       + ---      chr1     21-30       + |  0.511791
#>   [2]      chr1     11-20       + ---      chr1     51-55       + |  0.835552
#>   [3]      chr1     11-30       - ---      chr1     51-55       - |  0.708781
#>   [4]      chr1     11-30       - ---      chr2     51-60       - |  0.874206
#>              type
#>       <character>
#>   [1]         cis
#>   [2]         cis
#>   [3]         cis
#>   [4]       trans
#>   -------
#>   regions: 6 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_by("first") |> flank_left(-2) 
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors1
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1     11-12       + ---      chr1     21-30       + |  0.511791
#>   [2]      chr1     11-12       + ---      chr1     51-55       + |  0.835552
#>   [3]      chr1     11-12       - ---      chr1     51-55       - |  0.708781
#>   [4]      chr1     11-12       - ---      chr2     51-60       - |  0.874206
#>              type
#>       <character>
#>   [1]         cis
#>   [2]         cis
#>   [3]         cis
#>   [4]       trans
#>   -------
#>   regions: 6 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_by("second") |> flank_downstream(4)
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors2
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1     11-20       + ---      chr1     31-34       + |  0.511791
#>   [2]      chr1     11-20       + ---      chr1     56-59       + |  0.835552
#>   [3]      chr1     11-30       - ---      chr1     47-50       - |  0.708781
#>   [4]      chr1     11-30       - ---      chr2     47-50       - |  0.874206
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
# 2. Chained flanking of each set of anchors
####################################################################

gi |> 
  pin_by("first") |> flank_left(2) |> 
  pin_by("second") |> flank_right(2)
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors2
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1      9-10       + ---      chr1     31-32       + |  0.511791
#>   [2]      chr1      9-10       + ---      chr1     56-57       + |  0.835552
#>   [3]      chr1      9-10       - ---      chr1     56-57       - |  0.708781
#>   [4]      chr1      9-10       - ---      chr2     61-62       - |  0.874206
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
