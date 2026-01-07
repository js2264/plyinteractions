# Slice a GInteractions rows by their index

Slice a GInteractions rows by their index

## Usage

``` r
# S3 method for class 'GInteractions'
slice(.data, ...)
```

## Arguments

- .data:

  a GInteractions object

- ...:

  Integer indicating rows to keep.

## Value

a GInteractions object.

## Examples

``` r
gi <- read.table(text = "
chr1 1 10 chr1 1 10
chr2 1 10 chr2 1 10
chr3 1 10 chr3 1 10
chr4 1 10 chr4 1 10
chr5 1 10 chr5 1 10",
col.names = c(
    "seqnames1", "start1", "end1", 
    "seqnames2", "start2", "end2")
) |> 
  as_ginteractions()
  
####################################################################
# 1. Slice a GInteractions
####################################################################

gi |> slice(1, 2, 3)
#> GInteractions object with 3 interactions and 0 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle>
#>   [1]      chr1      1-10       * ---      chr1      1-10       *
#>   [2]      chr2      1-10       * ---      chr2      1-10       *
#>   [3]      chr3      1-10       * ---      chr3      1-10       *
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths
gi |> slice(-3)
#> GInteractions object with 4 interactions and 0 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle>
#>   [1]      chr1      1-10       * ---      chr1      1-10       *
#>   [2]      chr2      1-10       * ---      chr2      1-10       *
#>   [3]      chr4      1-10       * ---      chr4      1-10       *
#>   [4]      chr5      1-10       * ---      chr5      1-10       *
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths
gi |> slice(1:2, 5:4)
#> GInteractions object with 4 interactions and 0 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle>
#>   [1]      chr1      1-10       * ---      chr1      1-10       *
#>   [2]      chr2      1-10       * ---      chr2      1-10       *
#>   [3]      chr5      1-10       * ---      chr5      1-10       *
#>   [4]      chr4      1-10       * ---      chr4      1-10       *
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths
```
