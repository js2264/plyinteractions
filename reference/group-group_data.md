# GInteractions grouping metadata

GInteractions grouping metadata

## Usage

``` r
# S3 method for class 'GroupedGInteractions'
group_data(.data)

# S3 method for class 'GroupedGInteractions'
group_keys(.tbl, ...)

# S3 method for class 'GroupedGInteractions'
group_indices(.data, ...)

# S3 method for class 'GInteractions'
group_vars(x)

# S3 method for class 'GroupedGInteractions'
group_vars(x)

# S3 method for class 'GroupedGInteractions'
groups(x)

# S3 method for class 'GroupedGInteractions'
group_size(x)

# S3 method for class 'GroupedGInteractions'
n_groups(x)
```

## Arguments

- .data, .tbl, x:

  a GInteractions object

- ...:

  Ignored.

## Value

a GInteractions object.

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

ggi <- gi |> group_by(end1)
ggi
#> GroupedGInteractions object with 4 interactions and 2 metadata columns:
#> Groups: end1 [2]
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     11-20       * ---      chr1     21-30       * |         cis
#>   [2]      chr1     11-20       * ---      chr1     51-55       * |         cis
#>   [3]      chr1     11-30       * ---      chr1     51-55       * |         cis
#>   [4]      chr1     11-30       * ---      chr2     51-60       * |       trans
#>           score
#>       <numeric>
#>   [1]  0.805680
#>   [2]  0.814051
#>   [3]  0.403911
#>   [4]  0.218431
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
group_data(ggi)
#> DataFrame with 2 rows and 2 columns
#>        end1         .rows
#>   <integer> <IntegerList>
#> 1        20           1,2
#> 2        30           3,4
group_keys(ggi)
#> DataFrame with 2 rows and 1 column
#>        end1
#>   <integer>
#> 1        20
#> 2        30
group_rows(ggi)
#> IntegerList of length 2
#> [[1]] 1 2
#> [[2]] 3 4
group_indices(ggi)
#> integer-Rle of length 4 with 2 runs
#>   Lengths: 2 2
#>   Values : 1 2
group_vars(ggi)
#> [1] "end1"
groups(ggi)
#> [[1]]
#> end1
#> 
group_size(ggi)
#> [1] 2 2
n_groups(ggi)
#> [1] 2
```
