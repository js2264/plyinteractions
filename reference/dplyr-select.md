# Select columns within GInteractions metadata columns

Select columns within GInteractions metadata columns

## Usage

``` r
# S3 method for class 'GInteractions'
select(.data, ..., .drop_ranges = FALSE)
```

## Arguments

- .data:

  a GInteractions object

- ...:

  Integer indicating rows to keep.

- .drop_ranges:

  if TRUE, returns a DataFrame object. In this case, it enables
  selection of any column including core GInteractions columns.

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
  as_ginteractions() |> 
  mutate(score = runif(5)*100, cis = TRUE, gc = runif(5))
  
####################################################################
# 1. Select metadata columns from GInteractions by index
####################################################################

gi |> select(2, 1)
#> GInteractions object with 5 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |       cis
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <logical>
#>   [1]      chr1      1-10       * ---      chr1      1-10       * |      TRUE
#>   [2]      chr2      1-10       * ---      chr2      1-10       * |      TRUE
#>   [3]      chr3      1-10       * ---      chr3      1-10       * |      TRUE
#>   [4]      chr4      1-10       * ---      chr4      1-10       * |      TRUE
#>   [5]      chr5      1-10       * ---      chr5      1-10       * |      TRUE
#>           score
#>       <numeric>
#>   [1]   17.4676
#>   [2]   53.1574
#>   [3]   49.3637
#>   [4]   77.9309
#>   [5]   20.4178
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths
gi |> select(-3)
#> GInteractions object with 5 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1      1-10       * ---      chr1      1-10       * |   17.4676
#>   [2]      chr2      1-10       * ---      chr2      1-10       * |   53.1574
#>   [3]      chr3      1-10       * ---      chr3      1-10       * |   49.3637
#>   [4]      chr4      1-10       * ---      chr4      1-10       * |   77.9309
#>   [5]      chr5      1-10       * ---      chr5      1-10       * |   20.4178
#>             cis
#>       <logical>
#>   [1]      TRUE
#>   [2]      TRUE
#>   [3]      TRUE
#>   [4]      TRUE
#>   [5]      TRUE
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths

####################################################################
# 2. Select metadata columns from GInteractions by name
####################################################################

gi |> select(gc, score)
#> GInteractions object with 5 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        gc
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1      1-10       * ---      chr1      1-10       * | 0.7133973
#>   [2]      chr2      1-10       * ---      chr2      1-10       * | 0.0652161
#>   [3]      chr3      1-10       * ---      chr3      1-10       * | 0.3542068
#>   [4]      chr4      1-10       * ---      chr4      1-10       * | 0.8251994
#>   [5]      chr5      1-10       * ---      chr5      1-10       * | 0.2738182
#>           score
#>       <numeric>
#>   [1]   17.4676
#>   [2]   53.1574
#>   [3]   49.3637
#>   [4]   77.9309
#>   [5]   20.4178
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths

####################################################################
# 3. Select metadata columns from GInteractions with <tidy-select>
####################################################################

gi |> select(contains('s'))
#> GInteractions object with 5 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1      1-10       * ---      chr1      1-10       * |   17.4676
#>   [2]      chr2      1-10       * ---      chr2      1-10       * |   53.1574
#>   [3]      chr3      1-10       * ---      chr3      1-10       * |   49.3637
#>   [4]      chr4      1-10       * ---      chr4      1-10       * |   77.9309
#>   [5]      chr5      1-10       * ---      chr5      1-10       * |   20.4178
#>             cis
#>       <logical>
#>   [1]      TRUE
#>   [2]      TRUE
#>   [3]      TRUE
#>   [4]      TRUE
#>   [5]      TRUE
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths
gi |> select(matches('^s'))
#> GInteractions object with 5 interactions and 1 metadata column:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1      1-10       * ---      chr1      1-10       * |   17.4676
#>   [2]      chr2      1-10       * ---      chr2      1-10       * |   53.1574
#>   [3]      chr3      1-10       * ---      chr3      1-10       * |   49.3637
#>   [4]      chr4      1-10       * ---      chr4      1-10       * |   77.9309
#>   [5]      chr5      1-10       * ---      chr5      1-10       * |   20.4178
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths

####################################################################
# 4. Select core and metadata columns with .drop_ranges = TRUE
####################################################################

gi |> select(matches('^s'), .drop_ranges = TRUE)
#> DataFrame with 5 rows and 7 columns
#>   seqnames1    start1 strand1 seqnames2    start2 strand2     score
#>       <Rle> <integer>   <Rle>     <Rle> <integer>   <Rle> <numeric>
#> 1      chr1         1       *      chr1         1       *   17.4676
#> 2      chr2         1       *      chr2         1       *   53.1574
#> 3      chr3         1       *      chr3         1       *   49.3637
#> 4      chr4         1       *      chr4         1       *   77.9309
#> 5      chr5         1       *      chr5         1       *   20.4178
```
