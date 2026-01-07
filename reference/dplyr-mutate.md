# Mutate columns from a GInteractions object

Mutate columns from a GInteractions object

## Usage

``` r
# S3 method for class 'GInteractions'
mutate(.data, ...)
```

## Arguments

- .data:

  a GInteractions object

- ...:

  Optional named arguments specifying which the columns in .data to
  create/modify.

## Value

a GInteractions object.

## Examples

``` r
gi <- read.table(text = "
chr1 10 20 chr1 50 51
chr1 10 50 chr2 30 40",
col.names = c("chr1", "start1", "end1", "chr2", "start2", "end2")) |> 
  as_ginteractions(seqnames1 = chr1, seqnames2 = chr2)
  
####################################################################
# 1. Add metadata columns to a GInteractions object
####################################################################

gi |> 
  mutate(type = c('cis', 'trans'), score = runif(2)) |> 
  mutate(type2 = type)
#> GInteractions object with 2 interactions and 3 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     10-20       * ---      chr1     50-51       * |         cis
#>   [2]      chr1     10-50       * ---      chr2     30-40       * |       trans
#>           score       type2
#>       <numeric> <character>
#>   [1]  0.680163         cis
#>   [2]  0.498846       trans
#>   -------
#>   regions: 4 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

####################################################################
# 2. More complex, nested or inplace changes
####################################################################

gi |> 
  mutate(type = c('cis', 'trans'), score = runif(2)) |> 
  mutate(type2 = type) |> 
  mutate(count = c(1, 2), score = count * 2, new_col = paste0(type2, score))
#> GInteractions object with 2 interactions and 5 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     10-20       * ---      chr1     50-51       * |         cis
#>   [2]      chr1     10-50       * ---      chr2     30-40       * |       trans
#>           score       type2     count     new_col
#>       <numeric> <character> <numeric> <character>
#>   [1]         2         cis         1        cis2
#>   [2]         4       trans         2      trans4
#>   -------
#>   regions: 4 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

####################################################################
# 3. Core GInteractions columns can also be modified
####################################################################

gi |> 
  mutate(start1 = 1, end1 = 10, width2 = 30, strand2 = c('-', '+'))
#> GInteractions object with 2 interactions and 0 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle>
#>   [1]      chr1      1-10       * ---      chr1     50-79       -
#>   [2]      chr1      1-10       * ---      chr2     30-59       +
#>   -------
#>   regions: 3 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

# Note how the core columns are modified sequentially 

gi |> 
  mutate(start1 = 1, end1 = 10)
#> GInteractions object with 2 interactions and 0 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle>
#>   [1]      chr1      1-10       * ---      chr1     50-51       *
#>   [2]      chr1      1-10       * ---      chr2     30-40       *
#>   -------
#>   regions: 3 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> 
  mutate(start1 = 1, end1 = 10, width1 = 50)
#> GInteractions object with 2 interactions and 0 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle>
#>   [1]      chr1      1-50       * ---      chr1     50-51       *
#>   [2]      chr1      1-50       * ---      chr2     30-40       *
#>   -------
#>   regions: 3 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

####################################################################
# 4. Evaluating core GInteractions columns
####################################################################

gi |> 
  mutate(
    score = runif(2), 
    cis = seqnames1 == seqnames2, 
    distance = ifelse(cis, start2 - end1, NA)
  )
#> GInteractions object with 2 interactions and 3 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1     10-20       * ---      chr1     50-51       * | 0.0960242
#>   [2]      chr1     10-50       * ---      chr2     30-40       * | 0.7656002
#>         cis  distance
#>       <Rle> <integer>
#>   [1]  TRUE        30
#>   [2] FALSE      <NA>
#>   -------
#>   regions: 4 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
