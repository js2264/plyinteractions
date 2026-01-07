# Rename columns from a GInteractions with tidyverse-like `rename`

Rename columns from a GInteractions with tidyverse-like `rename`

## Usage

``` r
# S3 method for class 'GInteractions'
rename(.data, ...)
```

## Arguments

- .data:

  a GInteractions object

- ...:

  Use `new_name = old_name` to rename selected variables.

## Value

a GInteractions object.

## Examples

``` r
gi <- read.table(text = "
chr1 10 20 chr1 50 51
chr1 10 50 chr2 30 40",
col.names = c("chr1", "start1", "end1", "chr2", "start2", "end2")) |> 
  as_ginteractions(seqnames1 = chr1, seqnames2 = chr2) |> 
  mutate(type = c('cis', 'trans'), score = runif(2))
  
####################################################################
# 1. Rename metadata columns to a GInteractions object
####################################################################

gi |> rename(interaction_type = type, GC = score)
#> GInteractions object with 2 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> |
#>   [1]      chr1     10-20       * ---      chr1     50-51       * |
#>   [2]      chr1     10-50       * ---      chr2     30-40       * |
#>       interaction_type        GC
#>            <character> <numeric>
#>   [1]              cis  0.970521
#>   [2]            trans  0.389183
#>   -------
#>   regions: 4 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
