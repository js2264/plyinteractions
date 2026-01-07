# Group GInteractions by columns

Group GInteractions by columns

## Usage

``` r
# S3 method for class 'GInteractions'
group_by(.data, ..., .add = FALSE)

# S3 method for class 'DelegatingGInteractions'
group_by(.data, ..., .add = FALSE)

# S3 method for class 'GroupedGInteractions'
ungroup(x, ...)
```

## Arguments

- .data, x:

  a (Grouped)GInteractions object

- ...:

  Column(s) to group by.

- .add:

  When FALSE, the default, group_by() will override existing groups. To
  add to the existing groups, use .add = TRUE.

## Value

a `GroupedGInteractions` object. When a `(Anchored)PinnedGInteractions`
object is grouped, both anchoring and pinning are dropped.

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

####################################################################
# 1. Group by core column
####################################################################

gi |> group_by(end1)
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
#>   [1]  0.479025
#>   [2]  0.432171
#>   [3]  0.706434
#>   [4]  0.948577
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> group_by(end1, end2) |> group_data()
#> DataFrame with 4 rows and 3 columns
#>        end1      end2         .rows
#>   <integer> <integer> <IntegerList>
#> 1        20        30             1
#> 2        20        55             2
#> 3        30        55             3
#> 4        30        60             4

####################################################################
# 2. Group by metadata column
####################################################################

gi |> group_by(type) |> group_data()
#> DataFrame with 2 rows and 2 columns
#>          type         .rows
#>   <character> <IntegerList>
#> 1         cis         1,2,3
#> 2       trans             4

####################################################################
# 3. Combine core and metadata column grouping
####################################################################

gi |> group_by(end1, type)
#> GroupedGInteractions object with 4 interactions and 2 metadata columns:
#> Groups: end1, type [3]
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     11-20       * ---      chr1     21-30       * |         cis
#>   [2]      chr1     11-20       * ---      chr1     51-55       * |         cis
#>   [3]      chr1     11-30       * ---      chr1     51-55       * |         cis
#>   [4]      chr1     11-30       * ---      chr2     51-60       * |       trans
#>           score
#>       <numeric>
#>   [1]  0.479025
#>   [2]  0.432171
#>   [3]  0.706434
#>   [4]  0.948577
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
gi |> group_by(end1, type) |> group_data()
#> DataFrame with 3 rows and 3 columns
#>        end1        type         .rows
#>   <integer> <character> <IntegerList>
#> 1        20         cis           1,2
#> 2        30         cis             3
#> 3        30       trans             4

####################################################################
# 4. Create a new column and group by this new variable
####################################################################

gi |> group_by(class = c(1, 2, 1, 2))
#> GroupedGInteractions object with 4 interactions and 3 metadata columns:
#> Groups: class [2]
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     11-20       * ---      chr1     21-30       * |         cis
#>   [2]      chr1     11-20       * ---      chr1     51-55       * |         cis
#>   [3]      chr1     11-30       * ---      chr1     51-55       * |         cis
#>   [4]      chr1     11-30       * ---      chr2     51-60       * |       trans
#>           score     class
#>       <numeric> <numeric>
#>   [1]  0.479025         1
#>   [2]  0.432171         2
#>   [3]  0.706434         1
#>   [4]  0.948577         2
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

####################################################################
# 5. Replace or add groups to a GroupedGInteractions
####################################################################

ggi <- gi |> group_by(class = c(1, 2, 1, 2))
ggi |> group_data()
#> DataFrame with 2 rows and 2 columns
#>       class         .rows
#>   <numeric> <IntegerList>
#> 1         1           1,3
#> 2         2           2,4
ggi |> group_by(type) |> group_data()
#> DataFrame with 2 rows and 2 columns
#>          type         .rows
#>   <character> <IntegerList>
#> 1         cis         1,2,3
#> 2       trans             4
ggi |> group_by(type, .add = TRUE) |> group_data()
#> DataFrame with 3 rows and 3 columns
#>       class        type         .rows
#>   <numeric> <character> <IntegerList>
#> 1         1         cis           1,3
#> 2         2         cis             2
#> 3         2       trans             4

####################################################################
# 6. Ungroup GInteractions
####################################################################

ggi <- gi |> group_by(type, class = c(1, 2, 1, 2))
ggi
#> GroupedGInteractions object with 4 interactions and 3 metadata columns:
#> Groups: type, class [3]
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     11-20       * ---      chr1     21-30       * |         cis
#>   [2]      chr1     11-20       * ---      chr1     51-55       * |         cis
#>   [3]      chr1     11-30       * ---      chr1     51-55       * |         cis
#>   [4]      chr1     11-30       * ---      chr2     51-60       * |       trans
#>           score     class
#>       <numeric> <numeric>
#>   [1]  0.479025         1
#>   [2]  0.432171         2
#>   [3]  0.706434         1
#>   [4]  0.948577         2
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
ungroup(ggi, type)
#> GroupedGInteractions object with 4 interactions and 3 metadata columns:
#> Groups: class [2]
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     11-20       * ---      chr1     21-30       * |         cis
#>   [2]      chr1     11-20       * ---      chr1     51-55       * |         cis
#>   [3]      chr1     11-30       * ---      chr1     51-55       * |         cis
#>   [4]      chr1     11-30       * ---      chr2     51-60       * |       trans
#>           score     class
#>       <numeric> <numeric>
#>   [1]  0.479025         1
#>   [2]  0.432171         2
#>   [3]  0.706434         1
#>   [4]  0.948577         2
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
ungroup(ggi, class)
#> GroupedGInteractions object with 4 interactions and 3 metadata columns:
#> Groups: type [2]
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     11-20       * ---      chr1     21-30       * |         cis
#>   [2]      chr1     11-20       * ---      chr1     51-55       * |         cis
#>   [3]      chr1     11-30       * ---      chr1     51-55       * |         cis
#>   [4]      chr1     11-30       * ---      chr2     51-60       * |       trans
#>           score     class
#>       <numeric> <numeric>
#>   [1]  0.479025         1
#>   [2]  0.432171         2
#>   [3]  0.706434         1
#>   [4]  0.948577         2
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
