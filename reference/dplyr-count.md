# Count or tally GInteractions per group

Count or tally GInteractions per group

## Usage

``` r
# S3 method for class 'GroupedGInteractions'
tally(x, wt = NULL, sort = FALSE, name = NULL)

# S3 method for class 'GroupedGInteractions'
count(x, ..., wt = NULL, sort = FALSE, name = NULL)

# S3 method for class 'GInteractions'
count(x, ..., wt = NULL, sort = FALSE, name = NULL)
```

## Arguments

- x:

  A grouped GInteractions object

- wt:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Frequency weights. Can be `NULL` or a variable:

  - If `NULL` (the default), counts the number of rows in each group.

  - If a variable, computes `sum(wt)` for each group.

- sort:

  If `TRUE`, will show the largest groups at the top.

- name:

  The name of the new column in the output.

- ...:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Variables to group by.

## Value

a
`S4Vectors::`[`DataFrame()`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
object, with an added column with the count/tablly per group.

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
# 1. Tally groups
####################################################################

gi 
#> GInteractions object with 4 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1     11-20       + ---      chr1     21-30       + | 0.4035381
#>   [2]      chr1     11-20       + ---      chr1     51-55       + | 0.0636615
#>   [3]      chr1     11-30       - ---      chr1     51-55       - | 0.3887013
#>   [4]      chr1     11-30       - ---      chr2     51-60       - | 0.9755478
#>              type
#>       <character>
#>   [1]         cis
#>   [2]         cis
#>   [3]         cis
#>   [4]       trans
#>   -------
#>   regions: 6 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> group_by(strand1) |> tally()
#> DataFrame with 2 rows and 2 columns
#>   strand1         n
#>     <Rle> <integer>
#> 1       +         2
#> 2       -         2

gi |> group_by(type) |> tally()
#> DataFrame with 2 rows and 2 columns
#>          type         n
#>   <character> <integer>
#> 1         cis         3
#> 2       trans         1

gi |> group_by(type) |> tally(wt = score)
#> DataFrame with 2 rows and 2 columns
#>          type         n
#>   <character> <numeric>
#> 1         cis  0.855901
#> 2       trans  0.975548

####################################################################
# 2. Count per groups
####################################################################

gi |> count(type)
#> DataFrame with 2 rows and 2 columns
#>          type         n
#>   <character> <integer>
#> 1         cis         3
#> 2       trans         1

gi |> group_by(type) |> count(strand1)
#> DataFrame with 3 rows and 3 columns
#>          type strand1         n
#>   <character>   <Rle> <integer>
#> 1         cis       +         2
#> 2         cis       -         1
#> 3       trans       -         1

gi |> group_by(type, strand1) |> count(wt = score)
#> DataFrame with 3 rows and 3 columns
#>          type strand1         n
#>   <character>   <Rle> <numeric>
#> 1         cis       +  0.467200
#> 2         cis       -  0.388701
#> 3       trans       -  0.975548
```
