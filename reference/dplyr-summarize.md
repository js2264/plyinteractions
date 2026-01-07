# Summarize GInteractions per group

Summarize GInteractions per group

## Usage

``` r
# S3 method for class 'GroupedGInteractions'
summarise(.data, ...)

# S3 method for class 'GroupedGInteractions'
summarize(.data, ...)
```

## Arguments

- .data:

  a (grouped) GInteractions object

- ...:

  Name-value pairs of summary functions. The name will be the name of
  the variable in the result.

## Value

a
`S4Vectors::`[`DataFrame()`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
object:

- The rows come from the underlying `group_keys()`.

- The columns are a combination of the grouping keys and the summary
  expressions that you provide.

- GInteractions class is **not** preserved, as a call to `summarize`
  fundamentally creates a new data frame

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
# 1. Summarize a single column
####################################################################

gi
#> GInteractions object with 4 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1     11-20       + ---      chr1     21-30       + |  0.947764
#>   [2]      chr1     11-20       + ---      chr1     51-55       + |  0.542480
#>   [3]      chr1     11-30       - ---      chr1     51-55       - |  0.544603
#>   [4]      chr1     11-30       - ---      chr2     51-60       - |  0.278597
#>              type
#>       <character>
#>   [1]         cis
#>   [2]         cis
#>   [3]         cis
#>   [4]       trans
#>   -------
#>   regions: 6 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> group_by(type) |> summarize(m = mean(score))
#> DataFrame with 2 rows and 2 columns
#>          type         m
#>   <character> <numeric>
#> 1         cis  0.678283
#> 2       trans  0.278597

gi |> group_by(strand1) |> summarize(m = mean(score))
#> DataFrame with 2 rows and 2 columns
#>   strand1         m
#>     <Rle> <numeric>
#> 1       +  0.745122
#> 2       -  0.411600

df <- gi |> 
  group_by(strand1) |> 
  summarize(m = mean(score), n = table(seqnames2))
df
#> DataFrame with 2 rows and 3 columns
#>   strand1         m             n
#>     <Rle> <numeric> <IntegerList>
#> 1       +  0.745122           2,0
#> 2       -  0.411600           1,1

df$n
#> IntegerList of length 2
#> [["1"]] chr1=2 chr2=0
#> [["2"]] chr1=1 chr2=1

####################################################################
# 2. Summarize by multiple columns
####################################################################

gi |> 
  group_by(strand1, seqnames2) |> 
  summarise(m = mean(score), n = table(type))
#> DataFrame with 3 rows and 4 columns
#>   strand1 seqnames2         m             n
#>     <Rle>     <Rle> <numeric> <IntegerList>
#> 1       +      chr1  0.745122           2,0
#> 2       -      chr1  0.544603           1,0
#> 3       -      chr2  0.278597           0,1
```
