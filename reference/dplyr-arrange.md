# Arrange a GInteractions by a column

Arrange a GInteractions by a column

## Usage

``` r
# S3 method for class 'GInteractions'
arrange(.data, ...)
```

## Arguments

- .data:

  a GInteractions object

- ...:

  Variables, or functions of variables. Use dplyr::desc() to sort a
  variable in descending order.

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
#>   [1]  0.739944 0.7725215
#>   [2] 46.639350 0.8746007
#>   [3] 49.777739 0.1749406
#>   [4] 28.976724 0.0342413
#>   [5] 73.288199 0.3203857
#>   -------
#>   regions: 9 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths

####################################################################
# 1. Arrange GInteractions by a numerical column
####################################################################

gi |> arrange(gc)
#> GInteractions object with 5 interactions and 3 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |   cis
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <Rle>
#>   [1]      chr4      4-10       * ---      chr4      1-10       * |  TRUE
#>   [2]      chr3      3-10       * ---      chr3      1-10       * |  TRUE
#>   [3]      chr5      5-10       * ---      chr5      1-10       * |  TRUE
#>   [4]      chr1      1-10       * ---      chr1      1-10       * |  TRUE
#>   [5]      chr1      2-10       * ---      chr2      1-10       * | FALSE
#>           score        gc
#>       <numeric> <numeric>
#>   [1] 28.976724 0.0342413
#>   [2] 49.777739 0.1749406
#>   [3] 73.288199 0.3203857
#>   [4]  0.739944 0.7725215
#>   [5] 46.639350 0.8746007
#>   -------
#>   regions: 9 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths

####################################################################
# 2. Arrange GInteractions by a logical column
####################################################################

gi |> arrange(cis)
#> GInteractions object with 5 interactions and 3 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |   cis
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <Rle>
#>   [1]      chr1      2-10       * ---      chr2      1-10       * | FALSE
#>   [2]      chr1      1-10       * ---      chr1      1-10       * |  TRUE
#>   [3]      chr3      3-10       * ---      chr3      1-10       * |  TRUE
#>   [4]      chr4      4-10       * ---      chr4      1-10       * |  TRUE
#>   [5]      chr5      5-10       * ---      chr5      1-10       * |  TRUE
#>           score        gc
#>       <numeric> <numeric>
#>   [1] 46.639350 0.8746007
#>   [2]  0.739944 0.7725215
#>   [3] 49.777739 0.1749406
#>   [4] 28.976724 0.0342413
#>   [5] 73.288199 0.3203857
#>   -------
#>   regions: 9 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths

####################################################################
# 3. Arrange GInteractions by a factor
####################################################################

gi |> 
  mutate(rep = factor(c("rep1", "rep2", "rep1", "rep2", "rep1"))) |> 
  arrange(rep)
#> GInteractions object with 5 interactions and 4 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |   cis
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <Rle>
#>   [1]      chr1      1-10       * ---      chr1      1-10       * |  TRUE
#>   [2]      chr3      3-10       * ---      chr3      1-10       * |  TRUE
#>   [3]      chr5      5-10       * ---      chr5      1-10       * |  TRUE
#>   [4]      chr1      2-10       * ---      chr2      1-10       * | FALSE
#>   [5]      chr4      4-10       * ---      chr4      1-10       * |  TRUE
#>           score        gc      rep
#>       <numeric> <numeric> <factor>
#>   [1]  0.739944 0.7725215     rep1
#>   [2] 49.777739 0.1749406     rep1
#>   [3] 73.288199 0.3203857     rep1
#>   [4] 46.639350 0.8746007     rep2
#>   [5] 28.976724 0.0342413     rep2
#>   -------
#>   regions: 9 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths

####################################################################
# 4. Combine sorting variables
####################################################################

gi |> 
  mutate(rep = factor(c("rep1", "rep2", "rep1", "rep2", "rep1"))) |> 
  arrange(dplyr::desc(rep), score)
#> GInteractions object with 5 interactions and 4 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |   cis
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <Rle>
#>   [1]      chr4      4-10       * ---      chr4      1-10       * |  TRUE
#>   [2]      chr1      2-10       * ---      chr2      1-10       * | FALSE
#>   [3]      chr1      1-10       * ---      chr1      1-10       * |  TRUE
#>   [4]      chr3      3-10       * ---      chr3      1-10       * |  TRUE
#>   [5]      chr5      5-10       * ---      chr5      1-10       * |  TRUE
#>           score        gc      rep
#>       <numeric> <numeric> <factor>
#>   [1] 28.976724 0.0342413     rep2
#>   [2] 46.639350 0.8746007     rep2
#>   [3]  0.739944 0.7725215     rep1
#>   [4] 49.777739 0.1749406     rep1
#>   [5] 73.288199 0.3203857     rep1
#>   -------
#>   regions: 9 ranges and 0 metadata columns
#>   seqinfo: 5 sequences from an unspecified genome; no seqlengths
```
