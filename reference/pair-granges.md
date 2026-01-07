# Pairwise combination of a GRanges object

Create a GInteractions object from a GRanges object, containing all
possible entry pairs

## Usage

``` r
pair_granges(x)
```

## Arguments

- x:

  A GRanges object

## Value

A GInteractions object

## Examples

``` r
gr <- read.table(text = "
chr1 100 200
chr1 5000 5100
chr1 1000 5000
chr2 3000 3800",
col.names = c(
  "seqnames", "start", "end" 
)) |> plyranges::as_granges()

pair_granges(gr)
#> GInteractions object with 6 interactions and 0 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle>
#>   [1]      chr1   100-200       * ---      chr1 5000-5100       *
#>   [2]      chr1   100-200       * ---      chr1 1000-5000       *
#>   [3]      chr1   100-200       * ---      chr2 3000-3800       *
#>   [4]      chr1 5000-5100       * ---      chr1 1000-5000       *
#>   [5]      chr1 5000-5100       * ---      chr2 3000-3800       *
#>   [6]      chr1 1000-5000       * ---      chr2 3000-3800       *
#>   -------
#>   regions: 4 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
