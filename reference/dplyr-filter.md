# Subset a GInteractions with tidyverse-like `filter`

Subset a GInteractions with tidyverse-like `filter`

## Usage

``` r
# S3 method for class 'GInteractions'
filter(.data, ...)
```

## Arguments

- .data:

  a GInteractions object

- ...:

  Expressions that return a logical value, and are defined in terms of
  the variables in .data. If multiple expressions are included, they are
  combined with the & operator. Only rows for which all conditions
  evaluate to TRUE are kept.

## Value

a GInteractions object.

## Examples

``` r
gi <- read.table(text = "
chr1 1 10 chr1 1 10
chr1 2 10 chr2 1 10
chr3 3 10 chr3 1 10
chr4 4 10 chr4 1 10
chr5 5 10 chr5 1 10",
col.names = c(
    "seqnames1", "start1", "end1", 
    "seqnames2", "start2", "end2")
) |> 
  as_ginteractions() |> 
  mutate(cis = seqnames1 == seqnames2, score = runif(5)*100, gc = runif(5))
gi
#> GInteractions object with 5 interactions and 3 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |   cis
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <Rle>
#>   [1]      chr1      1-10       * ---      chr1      1-10       * |  TRUE
#>   [2]      chr1      2-10       * ---      chr2      1-10       * | FALSE
#>   [3]      chr3      3-10       * ---      chr3      1-10       * |  TRUE
#>   [4]      chr4      4-10       * ---      chr4      1-10       * |  TRUE
#>   [5]      chr5      5-10       * ---      chr5      1-10       * |  TRUE
#>           score        gc
#>       <numeric> <numeric>
#>   [1]  73.53196 0.5302125
#>   [2]  19.59567 0.6958239
#>   [3]  98.05397 0.6885560
#>   [4]  74.15215 0.0312303
#>   [5]   5.14463 0.2255625
#>   -------
#>   regions: 9 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths

####################################################################
# 1. Filter metadata columns from GInteractions by condition
####################################################################

gi |> filter(gc > 0.1)
#> GInteractions object with 4 interactions and 3 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |   cis
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <Rle>
#>   [1]      chr1      1-10       * ---      chr1      1-10       * |  TRUE
#>   [2]      chr1      2-10       * ---      chr2      1-10       * | FALSE
#>   [3]      chr3      3-10       * ---      chr3      1-10       * |  TRUE
#>   [4]      chr5      5-10       * ---      chr5      1-10       * |  TRUE
#>           score        gc
#>       <numeric> <numeric>
#>   [1]  73.53196  0.530212
#>   [2]  19.59567  0.695824
#>   [3]  98.05397  0.688556
#>   [4]   5.14463  0.225563
#>   -------
#>   regions: 9 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths
gi |> filter(gc > 0.1, score > 50)
#> GInteractions object with 2 interactions and 3 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |   cis
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <Rle>
#>   [1]      chr1      1-10       * ---      chr1      1-10       * |  TRUE
#>   [2]      chr3      3-10       * ---      chr3      1-10       * |  TRUE
#>           score        gc
#>       <numeric> <numeric>
#>   [1]    73.532  0.530212
#>   [2]    98.054  0.688556
#>   -------
#>   regions: 9 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths
gi |> filter(cis)
#> GInteractions object with 4 interactions and 3 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |   cis
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <Rle>
#>   [1]      chr1      1-10       * ---      chr1      1-10       * |  TRUE
#>   [2]      chr3      3-10       * ---      chr3      1-10       * |  TRUE
#>   [3]      chr4      4-10       * ---      chr4      1-10       * |  TRUE
#>   [4]      chr5      5-10       * ---      chr5      1-10       * |  TRUE
#>           score        gc
#>       <numeric> <numeric>
#>   [1]  73.53196 0.5302125
#>   [2]  98.05397 0.6885560
#>   [3]  74.15215 0.0312303
#>   [4]   5.14463 0.2255625
#>   -------
#>   regions: 9 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths

####################################################################
# 2. On-the-fly calculations
####################################################################

gi
#> GInteractions object with 5 interactions and 3 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |   cis
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <Rle>
#>   [1]      chr1      1-10       * ---      chr1      1-10       * |  TRUE
#>   [2]      chr1      2-10       * ---      chr2      1-10       * | FALSE
#>   [3]      chr3      3-10       * ---      chr3      1-10       * |  TRUE
#>   [4]      chr4      4-10       * ---      chr4      1-10       * |  TRUE
#>   [5]      chr5      5-10       * ---      chr5      1-10       * |  TRUE
#>           score        gc
#>       <numeric> <numeric>
#>   [1]  73.53196 0.5302125
#>   [2]  19.59567 0.6958239
#>   [3]  98.05397 0.6885560
#>   [4]  74.15215 0.0312303
#>   [5]   5.14463 0.2255625
#>   -------
#>   regions: 9 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths
gi |> filter(start1 >= start2 + 3)
#> GInteractions object with 2 interactions and 3 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |   cis
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <Rle>
#>   [1]      chr4      4-10       * ---      chr4      1-10       * |  TRUE
#>   [2]      chr5      5-10       * ---      chr5      1-10       * |  TRUE
#>           score        gc
#>       <numeric> <numeric>
#>   [1]  74.15215 0.0312303
#>   [2]   5.14463 0.2255625
#>   -------
#>   regions: 9 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths
gi |> filter(score * gc > score * 0.5)
#> GInteractions object with 3 interactions and 3 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |   cis
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <Rle>
#>   [1]      chr1      1-10       * ---      chr1      1-10       * |  TRUE
#>   [2]      chr1      2-10       * ---      chr2      1-10       * | FALSE
#>   [3]      chr3      3-10       * ---      chr3      1-10       * |  TRUE
#>           score        gc
#>       <numeric> <numeric>
#>   [1]   73.5320  0.530212
#>   [2]   19.5957  0.695824
#>   [3]   98.0540  0.688556
#>   -------
#>   regions: 9 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths
```
