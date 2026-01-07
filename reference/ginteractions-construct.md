# Construct a GInteractions object from a tibble, DataFrame or data.frame

The `as_ginteractions` function looks for column names in .data called
seqnames{1,2}, start{1,2}, end{1,2}, and strand{1,2} in order to
construct a GInteractions object. By default other columns in .data are
placed into the mcols (metadata columns) slot of the returned object.

## Usage

``` r
as_ginteractions(
  .data,
  ...,
  keep.extra.columns = TRUE,
  starts.in.df.are.0based = FALSE
)
```

## Arguments

- .data:

  A [`data.frame()`](https://rdrr.io/r/base/data.frame.html),
  [`S4Vectors::DataFrame()`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  or `tibble()` to construct a GInteractions object from.

- ...:

  Optional named arguments specifying which the columns in .data
  containin the core components a GInteractions object.

- keep.extra.columns:

  TRUE or FALSE (the default). If TRUE, the columns in df that are not
  used to form the genomic ranges of the returned GRanges object are
  then returned as metadata columns on the object. Otherwise, they are
  ignored.

- starts.in.df.are.0based:

  TRUE or FALSE (the default). If TRUE, then the start positions of the
  genomic ranges in df are considered to be 0-based and are converted to
  1-based in the returned GRanges object.

## Value

a GInteractions object.

## See also

[`InteractionSet::GInteractions()`](https://rdrr.io/pkg/InteractionSet/man/GInteractions-class.html)

## Examples

``` r
####################################################################
# 1. GInteractions from bedpe files imported into a data.frame
####################################################################

bedpe <- read.table(text = "
chr1 100 200 chr1 5000 5100 bedpe_example1 30 + -
chr1 1000 5000 chr1 3000 3800 bedpe_example2 100 + -",
col.names = c(
  "chrom1", "start1", "end1", 
  "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"))
bedpe |> 
  as_ginteractions(seqnames1 = chrom1, seqnames2 = chrom2)
#> GInteractions object with 2 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> |
#>   [1]      chr1   100-200       + ---      chr1 5000-5100       - |
#>   [2]      chr1 1000-5000       + ---      chr1 3000-3800       - |
#>                 name     score
#>          <character> <integer>
#>   [1] bedpe_example1        30
#>   [2] bedpe_example2       100
#>   -------
#>   regions: 4 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

####################################################################
# 2. GInteractions from standard pairs files imported into a data.frame
####################################################################

# Note how the pairs are 0-based and no "end" field is provided 
# (the standard pairs file format does not have "end" fields)
# We can provide width1 and width2 to fix this problem. 

pairs <- read.table(text = "
pair1 chr1 10000 chr1 20000 + +
pair2 chr1 50000 chr1 70000 + +
pair3 chr1 60000 chr2 10000 + +
pair4 chr1 30000 chr3 40000 + -", 
col.names = c(
  "pairID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2")
)
pairs |> 
  as_ginteractions(
    seqnames1 = chr1, start1 = pos1, width1 = 1000, 
    seqnames2 = chr2, start2 = pos2, width2 = 1000, 
    starts.in.df.are.0based = TRUE
  )
#> GInteractions object with 4 interactions and 1 metadata column:
#>       seqnames1     ranges1 strand1     seqnames2     ranges2 strand2 |
#>           <Rle>   <IRanges>   <Rle>         <Rle>   <IRanges>   <Rle> |
#>   [1]      chr1 10001-11000       + ---      chr1 20001-21000       + |
#>   [2]      chr1 50001-51000       + ---      chr1 70001-71000       + |
#>   [3]      chr1 60001-61000       + ---      chr2 10001-11000       + |
#>   [4]      chr1 30001-31000       + ---      chr3 40001-41000       - |
#>            pairID
#>       <character>
#>   [1]       pair1
#>   [2]       pair2
#>   [3]       pair3
#>   [4]       pair4
#>   -------
#>   regions: 8 ranges and 0 metadata columns
#>   seqinfo: 3 sequences from an unspecified genome; no seqlengths

####################################################################
# 3. GInteractions from data.frame with extra fields
####################################################################

df <- read.table(text = "
chr1 100 200 chr1 5000 5100
chr1 1000 5000 chr1 3000 3800",
col.names = c("chr1", "start1", "end1", "chr2", "start2", "end2"))
df |> 
  as_ginteractions(seqnames1 = chr1, seqnames2 = chr2)
#> GInteractions object with 2 interactions and 0 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle>
#>   [1]      chr1   100-200       * ---      chr1 5000-5100       *
#>   [2]      chr1 1000-5000       * ---      chr1 3000-3800       *
#>   -------
#>   regions: 4 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

df <- read.table(text = "
chr1 100 200 chr1 5000 5100
chr1 1000 5000 chr1 3000 3800",
col.names = c("chr1", "start1", "end1", "chr2", "start2", "end2"))
df |> 
  as_ginteractions(
    seqnames1 = chr1, seqnames2 = chr2, strand1 = '+', strand2 = '-'
  )
#> GInteractions object with 2 interactions and 0 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle>
#>   [1]      chr1   100-200       + ---      chr1 5000-5100       -
#>   [2]      chr1 1000-5000       + ---      chr1 3000-3800       -
#>   -------
#>   regions: 4 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

data.frame(type = "cis", count = 3) |> 
  as_ginteractions(
    seqnames1 = 'chr1', start1 = 1, end1 = 10,
    seqnames2 = 'chr1', start2 = 40, end2 = 50
  )
#> GInteractions object with 1 interaction and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1      1-10       * ---      chr1     40-50       * |         cis
#>           count
#>       <numeric>
#>   [1]         3
#>   -------
#>   regions: 2 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

####################################################################
# 4. GInteractions from a real like pairs files
####################################################################

pairsf <- system.file('extdata', 'pairs.gz', package = 'plyinteractions')
pairs <- read.table(pairsf, comment.char = '#', header = FALSE)
head(pairs)
#>                                           V1 V2  V3 V4     V5 V6 V7   V8   V9
#> 1  NS500150:527:HHGYNBGXF:3:21611:19085:3986 II 105 II  48548  +  - 1358 1681
#> 2  NS500150:527:HHGYNBGXF:4:13604:19734:2406 II 113 II  45003  -  + 1358 1658
#> 3 NS500150:527:HHGYNBGXF:2:11108:25178:11036 II 119 II 687251  -  + 1358 5550
#> 4   NS500150:527:HHGYNBGXF:1:22301:8468:1586 II 160 II  26124  +  - 1358 1510
#> 5  NS500150:527:HHGYNBGXF:4:23606:24037:2076 II 169 II  39052  +  + 1358 1613
#> 6  NS500150:527:HHGYNBGXF:1:12110:9220:19806 II 177 II  10285  +  - 1358 1416
pairs |> 
  as_ginteractions(
    seqnames1 = V2, start1 = V3, width1 = 1, strand1 = V6, 
    seqnames2 = V4, start2 = V5, width2 = 1, strand2 = V7,
    starts.in.df.are.0based = TRUE
  )
#> GInteractions object with 50000 interactions and 3 metadata columns:
#>           seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |
#>               <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> |
#>       [1]        II       106       + ---        II     48549       - |
#>       [2]        II       114       - ---        II     45004       + |
#>       [3]        II       120       - ---        II    687252       + |
#>       [4]        II       161       + ---        II     26125       - |
#>       [5]        II       170       + ---        II     39053       + |
#>       ...       ...       ...     ... ...       ...       ...     ... .
#>   [49996]        II     86997       + ---        II    487592       + |
#>   [49997]        II     86998       + ---        II     96354       - |
#>   [49998]        II     86998       + ---        II    114749       - |
#>   [49999]        II     86999       + ---        II     88956       + |
#>   [50000]        II     87000       + ---        II     87514       + |
#>                               V1        V8        V9
#>                      <character> <integer> <integer>
#>       [1] NS500150:527:HHGYNBG..      1358      1681
#>       [2] NS500150:527:HHGYNBG..      1358      1658
#>       [3] NS500150:527:HHGYNBG..      1358      5550
#>       [4] NS500150:527:HHGYNBG..      1358      1510
#>       [5] NS500150:527:HHGYNBG..      1358      1613
#>       ...                    ...       ...       ...
#>   [49996] NS500150:527:HHGYNBG..      1919      4378
#>   [49997] NS500150:527:HHGYNBG..      1919      1978
#>   [49998] NS500150:527:HHGYNBG..      1919      2092
#>   [49999] NS500150:527:HHGYNBG..      1919      1932
#>   [50000] NS500150:527:HHGYNBG..      1919      1922
#>   -------
#>   regions: 62911 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
