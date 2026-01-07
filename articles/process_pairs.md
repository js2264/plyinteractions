# Hi-C arithmetic with plyinteractions

The
*[plyinteractions](https://bioconductor.org/packages/3.23/plyinteractions)*
package facilitates data aggregation, for up to hundreds of thousands
and even millions of genomic interactions. In this vignette, we explore
several use cases which can arise when exploring Hi-C data stored in
`pairs` files.

We will use a real-life `pairs` file provided by the `4DN` Consortium.
This file has been generated from processing Hi-C performed in mouse
from brain cell primary culture during neural development (Bonev et al.,
Cell 2017). Pairs have been filtered to only those mapped over `chr13`.

``` r

library(plyinteractions)
#> Loading required package: InteractionSet
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff, setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind, colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique, unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse, colCounts, colCummaxs, colCummins, colCumprods, colCumsums, colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs, colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats, colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds, colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads, colWeightedMeans, colWeightedMedians, colWeightedSds, colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet, rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods, rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps, rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins, rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks, rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars, rowWeightedMads, rowWeightedMeans, rowWeightedMedians, rowWeightedSds, rowWeightedVars
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with 'browseVignettes()'. To cite Bioconductor, see 'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians
#> Loading required package: plyranges
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following object is masked from 'package:Biobase':
#> 
#>     combine
#> The following object is masked from 'package:matrixStats':
#> 
#>     count
#> The following objects are masked from 'package:GenomicRanges':
#> 
#>     intersect, setdiff, union
#> The following object is masked from 'package:Seqinfo':
#> 
#>     intersect
#> The following objects are masked from 'package:IRanges':
#> 
#>     collapse, desc, intersect, setdiff, slice, union
#> The following objects are masked from 'package:S4Vectors':
#> 
#>     first, intersect, rename, setdiff, setequal, union
#> The following objects are masked from 'package:BiocGenerics':
#> 
#>     combine, intersect, setdiff, setequal, union
#> The following object is masked from 'package:generics':
#> 
#>     explain
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> 
#> Attaching package: 'plyranges'
#> The following objects are masked from 'package:dplyr':
#> 
#>     between, n, n_distinct
#> 
#> Attaching package: 'plyinteractions'
#> The following objects are masked from 'package:plyranges':
#> 
#>     flank_downstream, flank_left, flank_right, flank_upstream, shift_downstream, shift_left, shift_right, shift_upstream
library(tidyverse)
#> ── Attaching core tidyverse packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 2.0.0 ──
#> ✔ forcats   1.0.1     ✔ readr     2.1.6
#> ✔ ggplot2   4.0.1     ✔ stringr   1.6.0
#> ✔ lubridate 1.9.4     ✔ tibble    3.3.0
#> ✔ purrr     1.2.0     ✔ tidyr     1.3.2
#> ── Conflicts ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ lubridate::%within%()   masks IRanges::%within%()
#> ✖ ggplot2::annotate()     masks plyinteractions::annotate()
#> ✖ plyranges::between()    masks dplyr::between()
#> ✖ dplyr::collapse()       masks IRanges::collapse()
#> ✖ dplyr::combine()        masks Biobase::combine(), BiocGenerics::combine()
#> ✖ dplyr::count()          masks matrixStats::count()
#> ✖ dplyr::desc()           masks IRanges::desc()
#> ✖ tidyr::expand()         masks S4Vectors::expand()
#> ✖ dplyr::filter()         masks stats::filter()
#> ✖ dplyr::first()          masks S4Vectors::first()
#> ✖ dplyr::lag()            masks stats::lag()
#> ✖ plyranges::n()          masks dplyr::n()
#> ✖ plyranges::n_distinct() masks dplyr::n_distinct()
#> ✖ ggplot2::Position()     masks BiocGenerics::Position(), base::Position()
#> ✖ purrr::reduce()         masks GenomicRanges::reduce(), IRanges::reduce()
#> ✖ dplyr::rename()         masks S4Vectors::rename()
#> ✖ lubridate::second()     masks S4Vectors::second()
#> ✖ lubridate::second<-()   masks S4Vectors::second<-()
#> ✖ dplyr::slice()          masks IRanges::slice()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

## Importing it in R
pairs_file <- HiContactsData::HiContactsData('mESCs', 'pairs.gz')
#> see ?HiContactsData and browseVignettes('HiContactsData') for documentation
#> loading from cache
pairs_df <- read.delim(
    pairs_file, sep = "\t", header = FALSE, comment.char = "#", nrows = 1e6
) |> 
    set_names(c(
        "ID", "seqnames1", "start1", 
        "seqnames2", "start2", "strand1", "strand2"
    )) 
pairs <- as_ginteractions(
    pairs_df, end1 = start1, end2 = start2, keep.extra.columns = TRUE
)
pairs
#> GInteractions object with 1000000 interactions and 1 metadata column:
#>             seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |                  ID
#>                 <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> |         <character>
#>         [1]     chr13  17057558       + ---     chr13  17176616       - |       SRR5339749.58
#>         [2]     chr13  68759440       - ---     chr13 113578864       - |      SRR5339749.105
#>         [3]     chr13  47940999       + ---     chr13  48134537       + |      SRR5339749.169
#>         [4]     chr13  80638451       + ---     chr13  80638826       - |      SRR5339749.170
#>         [5]     chr13   4362498       - ---     chr13  96982617       + |      SRR5339749.249
#>         ...       ...       ...     ... ...       ...       ...     ... .                 ...
#>    [999996]     chr13  17722638       - ---     chr13  20561010       - | SRR5339749.45907723
#>    [999997]     chr13  91988792       + ---     chr13  92333140       - | SRR5339749.45907730
#>    [999998]     chr13  98284294       - ---     chr13  98531510       + | SRR5339749.45907742
#>    [999999]     chr13  36053113       + ---     chr13  36076689       - | SRR5339749.45907855
#>   [1000000]     chr13  25585835       + ---     chr13  26197381       - | SRR5339749.45907876
#>   -------
#>   regions: 1936743 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

## Estimating pairs filtering thresholds

We can first *in silico* digest the mouse genome to obtain the
coordinates of each genomic fragment after digestion by **DpnII and
HinfI**.

``` r

## Prepare DpnII/HinfI-digested genomic fragments
library(Biostrings)
#> Loading required package: XVector
#> 
#> Attaching package: 'XVector'
#> The following object is masked from 'package:purrr':
#> 
#>     compact
#> 
#> Attaching package: 'Biostrings'
#> The following object is masked from 'package:base':
#> 
#>     strsplit
genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
cutter <- DNAStringSet(c("GATC", "GANTC"))  ## DpnII/HinfI cutting site
fragments <- BiocParallel::bplapply(BPPARAM = BiocParallel::MulticoreParam(workers = 2), 
    names(genome), function(.x) {
        seq <- genome[[.x]]
        mids <- lapply(
            cutter, 
            function(cutsite) {
                hits <- matchPattern(cutsite, seq, fixed = "subject")
                start(hits) + {end(hits) - start(hits)}
            }
        ) |> unlist() |> sort()
        GRanges(seqnames = .x, IRanges(
            start = c(1, mids), end = c(mids-1, length(seq))
        ))
    }
) |> 
    set_names(names(genome)) |> 
    GRangesList() |> 
    unlist()
fragments$binID <- seq_along(fragments)
```

We can then use the
[`annotate()`](../reference/ginteractions-annotate.md) function from
*[plyinteractions](https://bioconductor.org/packages/3.23/plyinteractions)*
to recover, for each interaction, which restriction enzyme fragment each
anchor overlaps with, and how many restriction enzyme cutting sites are
found between them.

``` r

## Annotate for each anchor set which genomic fragment it overlaps with
annotated_pairs <- pairs |> 
    plyinteractions::annotate(fragments, by = "binID") |> 
    mutate(n_fragments = binID.2 - binID.1, group = paste0(strand1, strand2))
annotated_pairs
#> GInteractions object with 1000000 interactions and 5 metadata columns:
#>             seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |                  ID   binID.1   binID.2 n_fragments       group
#>                 <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> |         <character> <integer> <integer>   <integer> <character>
#>         [1]     chr13  17057558       + ---     chr13  17176616       - |       SRR5339749.58   9591352   9592012         660          +-
#>         [2]     chr13  68759440       - ---     chr13 113578864       - |      SRR5339749.105   9880169  10124404      244235          --
#>         [3]     chr13  47940999       + ---     chr13  48134537       + |      SRR5339749.169   9762274   9763393        1119          ++
#>         [4]     chr13  80638451       + ---     chr13  80638826       - |      SRR5339749.170   9946878   9946878           0          +-
#>         [5]     chr13   4362498       - ---     chr13  96982617       + |      SRR5339749.249   9521271  10034142      512871          -+
#>         ...       ...       ...     ... ...       ...       ...     ... .                 ...       ...       ...         ...         ...
#>    [999996]     chr13  17722638       - ---     chr13  20561010       - | SRR5339749.45907723   9594967   9610226       15259          --
#>    [999997]     chr13  91988792       + ---     chr13  92333140       - | SRR5339749.45907730  10007449  10009172        1723          +-
#>    [999998]     chr13  98284294       - ---     chr13  98531510       + | SRR5339749.45907742  10041410  10042771        1361          -+
#>    [999999]     chr13  36053113       + ---     chr13  36076689       - | SRR5339749.45907855   9696787   9696912         125          +-
#>   [1000000]     chr13  25585835       + ---     chr13  26197381       - | SRR5339749.45907876   9638511   9641687        3176          +-
#>   -------
#>   regions: 1936743 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

Next, we can plot the distribution of `strand1` and `strand2`
cominations as a function of the number of restriction enzyme cutting
sites between anchors of each interaction.

``` r

df <- annotated_pairs |> 
    head(n = 1e6) |> 
    group_by(strand1, strand2, n_fragments) |> 
    count() |> 
    as_tibble() |> 
    mutate(group = paste0(strand1, strand2)) |> 
    select(group, n_fragments, n)
ggplot(df, aes(x = n_fragments, y = n, group = group, col = group)) + 
    geom_line() + 
    geom_point() + 
    xlim(c(0, 15)) + 
    annotation_logticks(sides = 'l') + 
    theme_bw() + 
    labs(
        x = "Number of restriction sites between anchors", 
        y = "Number of pairs"
    )
#> Warning: Removed 267493 rows containing missing values or values outside the scale range (`geom_line()`).
#> Warning: Removed 267493 rows containing missing values or values outside the scale range (`geom_point()`).
```

![](process_pairs_files/figure-html/getDistrib-1.png)

From this distribution, we can see that `--` and `++` pairs have a
decreasing frequency over increasing numbers of cut sites between
anchors of each interaction. These pairs are unambiguous, as the
orientation of each sequenced end can only come from true cutting and
religation event, (except the set of `--` and `++` pairs which have `0`
cut sites between each anchor, which cannot be explained); all these
pairs can be kept.

The over-representation of `+-` pairs at short distance likely represent
uncut fragments subsequently sequenced on each end. The
under-representation of `-+` pairs at short distance likely represent
self-religated fragments. We can estimate a threshold for each of these
pairs sets by computing the MAD and expected , as described in [Cournac
et al.,
2012](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-436).

``` r

filters <- df |> 
    filter(n_fragments <= 50) |> 
    arrange(n_fragments) |> 
    group_by(n_fragments) |> 
    mutate(median = median(n)) |> 
    ungroup() |> 
    mutate(MAD = median(abs(n - median))) |> 
    mutate(withinMAD = abs(n - median) <= MAD / 0.67449) |> 
    filter(withinMAD) |> 
    slice_head(by = group, n = 1) |> 
    select(group, n_fragments) |> 
    rename(threshold = n_fragments)
filters
#> # A tibble: 4 × 2
#>   group threshold
#>   <chr>     <int>
#> 1 ++            1
#> 2 --            1
#> 3 -+            8
#> 4 +-           10
```

## Filtering pairs using appropriate thresholds

``` r

annotated_pairs <- annotated_pairs |> 
    mutate(threshold = left_join(as_tibble(mcols(annotated_pairs)), filters)$threshold) |> 
    mutate(type = case_when(
        group %in% c('--', '++') & n_fragments < threshold ~ "excluded", 
        group == '+-' & n_fragments < threshold ~ "uncut", 
        group == '-+' & n_fragments < threshold ~ "religated", 
        .default = "kept"
    ))
#> Joining with `by = join_by(group)`
mcols(annotated_pairs) |>
    as_tibble() |> 
    count(type) |> 
    mutate(n = scales::percent(n/sum(n)))
#> # A tibble: 4 × 2
#>   type      n     
#>   <chr>     <chr> 
#> 1 excluded  1.08% 
#> 2 kept      78.70%
#> 3 religated 0.40% 
#> 4 uncut     19.82%

filtered_pairs <- filter(annotated_pairs, type == 'kept')
```

## Computing distance law from pairs

Another typical step when analyzing Hi-C processed data is the modeling
of a so-called “distance law”, (a.k.a “P(s)”), which describes the
genomic distance-dependent contact frequency between pairs of genomic
loci from a Hi-C experiment.

We can easily recover the distance between the two anchors of each
interaction (noted *s*) and plot the interaction frequency (noted
*P(s)*) as a function of this genomic distance.

### Plotting distance law: first try

``` r

dat <- filtered_pairs |> 
    mutate(s = abs(end2 - start1)) |> 
    group_by(s) |> 
    count(name = "n") |>
    as_tibble() |> 
    mutate(Ps = n/sum(n)) 
p <- ggplot(dat, aes(x = s, y = Ps)) + geom_line()
p
```

![](process_pairs_files/figure-html/Ps1-1.png)

This is not very informative, as the distances span several orders of
magnitude in both dimensions.

### Second try: switching to logarithmic scale

Switching to a `log` scale in
*[ggplot2](https://CRAN.R-project.org/package=ggplot2)* is very easy.

``` r

p + scale_x_log10() + scale_y_log10() + annotation_logticks()
```

![](process_pairs_files/figure-html/Ps2-1.png)

### Third try: aggregating data before plotting

The previous P(s) plot is precise at the base-pair resolution. We can
aggregate counts by binned distances:

``` r

# Calculate distance breaks evenly spaced on a log scale (base 1.1)
x <- 1.1^(1:200-1)
lmc <- coef(lm(c(1,1161443398)~c(x[1], x[200])))
bins_breaks <- unique(round(lmc[2]*x + lmc[1]))
bins_widths <- lead(bins_breaks) - bins_breaks

# Bin distances
dat <- filtered_pairs |> 
    mutate(s = abs(end2 - start1)) |> 
    mutate(
        binned_s = bins_breaks[as.numeric(cut(s, bins_breaks))], 
        bin_width = bins_widths[as.numeric(cut(s, bins_breaks))]
    ) |> 
    group_by(binned_s, bin_width) |> 
    count(name = "n") |>
    as_tibble() |> 
    mutate(Ps = n / sum(n) / bin_width)

# Plot results
ggplot(dat, aes(x = binned_s, y = Ps)) + geom_line() + 
    scale_x_log10() + scale_y_log10() + annotation_logticks()
```

![](process_pairs_files/figure-html/Ps3-1.png)

### With some polishing

``` r

ggplot(dat, aes(x = binned_s, y = Ps)) + 
    geom_line() + 
    scale_x_log10(limits = c(1e3, 1e8)) +    ## This changes x axis to log scale
    scale_y_log10() +                        ## This changes y axis to log scale
    annotation_logticks() +                  ## This adds log ticks
    labs(
        x = "Genomic distance (s)", 
        y = "P(s)", 
        title = "Distance-dependent genomic frequency P(s) in mESC (chr. 13)"
    ) +                                      ## This fixes axes titles
    theme_bw()                               ## This changes default plot theme
#> Warning: Removed 41 rows containing missing values or values outside the scale range (`geom_line()`).
```

![](process_pairs_files/figure-html/Ps4-1.png)

## Reproducibility

`R` session information:

    #> ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
    #>  setting  value
    #>  version  R Under development (unstable) (2026-01-03 r89269)
    #>  os       Ubuntu 24.04.3 LTS
    #>  system   x86_64, linux-gnu
    #>  ui       X11
    #>  language en
    #>  collate  en_US.UTF-8
    #>  ctype    en_US.UTF-8
    #>  tz       UTC
    #>  date     2026-01-07
    #>  pandoc   3.8.3 @ /usr/bin/ (via rmarkdown)
    #>  quarto   1.8.26 @ /usr/local/bin/quarto
    #> 
    #> ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
    #>  package                      * version   date (UTC) lib source
    #>  abind                          1.4-8     2024-09-12 [1] CRAN (R 4.6.0)
    #>  AnnotationDbi                  1.73.0    2025-10-31 [1] Bioconductor 3.23 (R 4.6.0)
    #>  AnnotationHub                * 4.1.0     2025-10-31 [1] Bioconductor 3.23 (R 4.6.0)
    #>  Biobase                      * 2.71.0    2025-10-30 [1] Bioconductor 3.23 (R 4.6.0)
    #>  BiocFileCache                * 3.1.0     2025-10-30 [1] Bioconductor 3.23 (R 4.6.0)
    #>  BiocGenerics                 * 0.57.0    2025-10-30 [1] Bioconductor 3.23 (R 4.6.0)
    #>  BiocIO                         1.21.0    2025-10-30 [1] Bioconductor 3.23 (R 4.6.0)
    #>  BiocManager                    1.30.27   2025-11-14 [1] CRAN (R 4.6.0)
    #>  BiocParallel                   1.45.0    2025-10-30 [1] Bioconductor 3.23 (R 4.6.0)
    #>  BiocStyle                    * 2.39.0    2025-10-30 [1] Bioconductor 3.23 (R 4.6.0)
    #>  BiocVersion                    3.23.1    2025-10-30 [2] Bioconductor 3.23 (R 4.6.0)
    #>  Biostrings                   * 2.79.3    2025-12-17 [1] Bioconductor 3.23 (R 4.6.0)
    #>  bit                            4.6.0     2025-03-06 [1] CRAN (R 4.6.0)
    #>  bit64                          4.6.0-1   2025-01-16 [1] CRAN (R 4.6.0)
    #>  bitops                         1.0-9     2024-10-03 [1] CRAN (R 4.6.0)
    #>  blob                           1.2.4     2023-03-17 [1] CRAN (R 4.6.0)
    #>  bookdown                       0.46      2025-12-05 [1] CRAN (R 4.6.0)
    #>  BSgenome                       1.79.1    2025-11-04 [1] Bioconductor 3.23 (R 4.6.0)
    #>  BSgenome.Mmusculus.UCSC.mm10   1.4.3     2025-12-12 [1] Bioconductor
    #>  bslib                          0.9.0     2025-01-30 [2] CRAN (R 4.6.0)
    #>  cachem                         1.1.0     2024-05-16 [2] CRAN (R 4.6.0)
    #>  cigarillo                      1.1.0     2025-10-31 [1] Bioconductor 3.23 (R 4.6.0)
    #>  cli                            3.6.5     2025-04-23 [2] CRAN (R 4.6.0)
    #>  codetools                      0.2-20    2024-03-31 [3] CRAN (R 4.6.0)
    #>  crayon                         1.5.3     2024-06-20 [2] CRAN (R 4.6.0)
    #>  curl                           7.0.0     2025-08-19 [2] CRAN (R 4.6.0)
    #>  DBI                            1.2.3     2024-06-02 [1] CRAN (R 4.6.0)
    #>  dbplyr                       * 2.5.1     2025-09-10 [1] CRAN (R 4.6.0)
    #>  DelayedArray                   0.37.0    2025-10-31 [1] Bioconductor 3.23 (R 4.6.0)
    #>  desc                           1.4.3     2023-12-10 [2] CRAN (R 4.6.0)
    #>  digest                         0.6.39    2025-11-19 [2] CRAN (R 4.6.0)
    #>  dplyr                        * 1.1.4     2023-11-17 [1] CRAN (R 4.6.0)
    #>  evaluate                       1.0.5     2025-08-27 [2] CRAN (R 4.6.0)
    #>  ExperimentHub                * 3.1.0     2025-10-31 [1] Bioconductor 3.23 (R 4.6.0)
    #>  farver                         2.1.2     2024-05-13 [1] CRAN (R 4.6.0)
    #>  fastmap                        1.2.0     2024-05-15 [2] CRAN (R 4.6.0)
    #>  filelock                       1.0.3     2023-12-11 [1] CRAN (R 4.6.0)
    #>  forcats                      * 1.0.1     2025-09-25 [1] CRAN (R 4.6.0)
    #>  fs                             1.6.6     2025-04-12 [2] CRAN (R 4.6.0)
    #>  generics                     * 0.1.4     2025-05-09 [1] CRAN (R 4.6.0)
    #>  GenomicAlignments              1.47.0    2025-10-31 [1] Bioconductor 3.23 (R 4.6.0)
    #>  GenomicRanges                * 1.63.1    2025-12-08 [1] Bioconductor 3.23 (R 4.6.0)
    #>  ggplot2                      * 4.0.1     2025-11-14 [1] CRAN (R 4.6.0)
    #>  glue                           1.8.0     2024-09-30 [2] CRAN (R 4.6.0)
    #>  gtable                         0.3.6     2024-10-25 [1] CRAN (R 4.6.0)
    #>  HiContactsData               * 1.13.0    2025-11-04 [1] Bioconductor 3.23 (R 4.6.0)
    #>  hms                            1.1.4     2025-10-17 [1] CRAN (R 4.6.0)
    #>  htmltools                      0.5.9     2025-12-04 [2] CRAN (R 4.6.0)
    #>  htmlwidgets                    1.6.4     2023-12-06 [2] CRAN (R 4.6.0)
    #>  httr                           1.4.7     2023-08-15 [1] CRAN (R 4.6.0)
    #>  httr2                          1.2.2     2025-12-08 [2] CRAN (R 4.6.0)
    #>  InteractionSet               * 1.39.0    2025-10-31 [1] Bioconductor 3.23 (R 4.6.0)
    #>  IRanges                      * 2.45.0    2025-10-31 [1] Bioconductor 3.23 (R 4.6.0)
    #>  jquerylib                      0.1.4     2021-04-26 [2] CRAN (R 4.6.0)
    #>  jsonlite                       2.0.0     2025-03-27 [2] CRAN (R 4.6.0)
    #>  KEGGREST                       1.51.1    2025-11-17 [1] Bioconductor 3.23 (R 4.6.0)
    #>  knitr                          1.51      2025-12-20 [2] CRAN (R 4.6.0)
    #>  labeling                       0.4.3     2023-08-29 [1] CRAN (R 4.6.0)
    #>  lattice                        0.22-7    2025-04-02 [3] CRAN (R 4.6.0)
    #>  lifecycle                      1.0.4     2023-11-07 [2] CRAN (R 4.6.0)
    #>  lubridate                    * 1.9.4     2024-12-08 [1] CRAN (R 4.6.0)
    #>  magrittr                       2.0.4     2025-09-12 [2] CRAN (R 4.6.0)
    #>  Matrix                         1.7-4     2025-08-28 [3] CRAN (R 4.6.0)
    #>  MatrixGenerics               * 1.23.0    2025-10-30 [1] Bioconductor 3.23 (R 4.6.0)
    #>  matrixStats                  * 1.5.0     2025-01-07 [1] CRAN (R 4.6.0)
    #>  memoise                        2.0.1     2021-11-26 [2] CRAN (R 4.6.0)
    #>  otel                           0.2.0     2025-08-29 [2] CRAN (R 4.6.0)
    #>  pillar                         1.11.1    2025-09-17 [2] CRAN (R 4.6.0)
    #>  pkgconfig                      2.0.3     2019-09-22 [2] CRAN (R 4.6.0)
    #>  pkgdown                        2.2.0     2025-11-06 [1] CRAN (R 4.6.0)
    #>  plyinteractions              * 1.9.1     2026-01-07 [1] Bioconductor
    #>  plyranges                    * 1.31.1    2025-11-07 [1] Bioconductor 3.23 (R 4.6.0)
    #>  png                            0.1-8     2022-11-29 [1] CRAN (R 4.6.0)
    #>  purrr                        * 1.2.0     2025-11-04 [2] CRAN (R 4.6.0)
    #>  R6                             2.6.1     2025-02-15 [2] CRAN (R 4.6.0)
    #>  ragg                           1.5.0     2025-09-02 [2] CRAN (R 4.6.0)
    #>  rappdirs                       0.3.3     2021-01-31 [2] CRAN (R 4.6.0)
    #>  RColorBrewer                   1.1-3     2022-04-03 [1] CRAN (R 4.6.0)
    #>  Rcpp                           1.1.0.8.1 2025-12-08 [2] CRAN (R 4.6.0)
    #>  RCurl                          1.98-1.17 2025-03-22 [1] CRAN (R 4.6.0)
    #>  readr                        * 2.1.6     2025-11-14 [1] CRAN (R 4.6.0)
    #>  restfulr                       0.0.16    2025-06-27 [1] CRAN (R 4.6.0)
    #>  rjson                          0.2.23    2024-09-16 [1] CRAN (R 4.6.0)
    #>  rlang                          1.1.6     2025-04-11 [2] CRAN (R 4.6.0)
    #>  rmarkdown                      2.30      2025-09-28 [1] CRAN (R 4.6.0)
    #>  Rsamtools                      2.27.0    2025-10-31 [1] Bioconductor 3.23 (R 4.6.0)
    #>  RSQLite                        2.4.5     2025-11-30 [1] CRAN (R 4.6.0)
    #>  rtracklayer                    1.71.3    2025-12-14 [1] Bioconductor 3.23 (R 4.6.0)
    #>  S4Arrays                       1.11.1    2025-11-25 [1] Bioconductor 3.23 (R 4.6.0)
    #>  S4Vectors                    * 0.49.0    2025-10-30 [1] Bioconductor 3.23 (R 4.6.0)
    #>  S7                             0.2.1     2025-11-14 [1] CRAN (R 4.6.0)
    #>  sass                           0.4.10    2025-04-11 [2] CRAN (R 4.6.0)
    #>  scales                         1.4.0     2025-04-24 [1] CRAN (R 4.6.0)
    #>  Seqinfo                      * 1.1.0     2025-10-31 [1] Bioconductor 3.23 (R 4.6.0)
    #>  sessioninfo                  * 1.2.3     2025-02-05 [2] CRAN (R 4.6.0)
    #>  SparseArray                    1.11.10   2025-12-16 [1] Bioconductor 3.23 (R 4.6.0)
    #>  stringi                        1.8.7     2025-03-27 [2] CRAN (R 4.6.0)
    #>  stringr                      * 1.6.0     2025-11-04 [2] CRAN (R 4.6.0)
    #>  SummarizedExperiment         * 1.41.0    2025-10-31 [1] Bioconductor 3.23 (R 4.6.0)
    #>  systemfonts                    1.3.1     2025-10-01 [2] CRAN (R 4.6.0)
    #>  textshaping                    1.0.4     2025-10-10 [2] CRAN (R 4.6.0)
    #>  tibble                       * 3.3.0     2025-06-08 [2] CRAN (R 4.6.0)
    #>  tidyr                        * 1.3.2     2025-12-19 [1] CRAN (R 4.6.0)
    #>  tidyselect                     1.2.1     2024-03-11 [1] CRAN (R 4.6.0)
    #>  tidyverse                    * 2.0.0     2023-02-22 [1] CRAN (R 4.6.0)
    #>  timechange                     0.3.0     2024-01-18 [1] CRAN (R 4.6.0)
    #>  tzdb                           0.5.0     2025-03-15 [1] CRAN (R 4.6.0)
    #>  utf8                           1.2.6     2025-06-08 [2] CRAN (R 4.6.0)
    #>  vctrs                          0.6.5     2023-12-01 [2] CRAN (R 4.6.0)
    #>  withr                          3.0.2     2024-10-28 [2] CRAN (R 4.6.0)
    #>  xfun                           0.55      2025-12-16 [2] CRAN (R 4.6.0)
    #>  XML                            3.99-0.20 2025-11-08 [1] CRAN (R 4.6.0)
    #>  XVector                      * 0.51.0    2025-10-31 [1] Bioconductor 3.23 (R 4.6.0)
    #>  yaml                           2.3.12    2025-12-10 [2] CRAN (R 4.6.0)
    #> 
    #>  [1] /__w/_temp/Library
    #>  [2] /usr/local/lib/R/site-library
    #>  [3] /usr/local/lib/R/library
    #>  * ── Packages attached to the search path.
    #> 
    #> ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
