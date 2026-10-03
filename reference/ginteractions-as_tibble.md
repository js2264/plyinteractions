# Turn a GInteractions object into a tibble

`as_tibble()` returns the interactions of a GInteractions object as a
tibble: one row per interaction, with the coordinates of both anchors
(`seqnames1`, `start1`, `end1`, `width1`, `strand1`, then the same for
the second anchor) followed by the metadata columns. It is the reverse
of [`as_ginteractions()`](ginteractions-construct.md).

## Usage

``` r
# S3 method for class 'GInteractions'
as_tibble(x, ...)
```

## Arguments

- x:

  A GInteractions object.

- ...:

  Passed to
  [`tibble::as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html).

## Value

A tibble.

## Examples

``` r
gi <- read.table(text = "
chr1 11 20 chr1 21 30 + +
chr1 11 20 chr1 51 55 + +
chr1 11 30 chr2 51 60 - -",
col.names = c(
  "seqnames1", "start1", "end1", 
  "seqnames2", "start2", "end2", "strand1", "strand2")
) |> 
  as_ginteractions()
gi$type <- c("cis", "cis", "trans")
as_tibble(gi)
#> # A tibble: 3 × 11
#>   seqnames1 start1  end1 width1 strand1 seqnames2 start2  end2 width2 strand2
#>   <fct>      <int> <int>  <int> <fct>   <fct>      <int> <int>  <int> <fct>  
#> 1 chr1          11    20     10 +       chr1          21    30     10 +      
#> 2 chr1          11    20     10 +       chr1          51    55      5 +      
#> 3 chr1          11    30     20 -       chr2          51    60     10 -      
#> # ℹ 1 more variable: type <chr>
```
