# show method for `GInteractions` objects

show method for `GInteractions` objects

## Arguments

- object:

  a `(Anchored/Pinned/Grouped)GInteractions` object

## Value

`Prints a message to the console describing the contents of a `GInteractions\`
object.

## Examples

``` r
pairsf <- system.file('extdata', 'pairs.gz', package = 'plyinteractions')
pairs <- read.table(pairsf, comment.char = '#', header = FALSE)
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
