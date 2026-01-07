# Annotate both anchors of a GInteractions

For each interaction in a `GInteractions` object, `annotate` returns the
pairs of annotations from the `GRanges` object it overlaps with.

## Usage

``` r
annotate(x, y, by)

annotate_directed(x, y, by)

# S4 method for class 'GInteractions,GRanges,character'
annotate(x, y, by)

# S4 method for class 'GInteractions,GRanges,character'
annotate_directed(x, y, by)
```

## Arguments

- x:

  a GInteractions object

- y:

  a GRanges object to extract annotations from

- by:

  Column name from `y` to use to extract annotations

## Value

a GInteractions object with two extra metadata columns named `by.1` and
`by.2`.

## Examples

``` r
####################################################################
# 1. Basic example
####################################################################

gi <- read.table(text = "  
    chr1 11 20 - chr1 21 30 + 
    chr1 21 30 + chr2 51 60 +",  
    col.names = c(
        "seqnames1", "start1", "end1", "strand1", 
        "seqnames2", "start2", "end2", "strand2"
    )
) |> as_ginteractions() 

gr <- GenomicRanges::GRanges(c("chr1:20-30:+", "chr2:55-65:+")) |>
    mutate(id = 1:2)

annotate(gi, gr, by = 'id')
#> GInteractions object with 2 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |      id.1
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>   [1]      chr1     11-20       - ---      chr1     21-30       + |         1
#>   [2]      chr1     21-30       + ---      chr2     51-60       + |         1
#>            id.2
#>       <integer>
#>   [1]         1
#>   [2]         2
#>   -------
#>   regions: 3 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

annotate_directed(gi, gr, by = 'id')
#> GInteractions object with 2 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |      id.1
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <integer>
#>   [1]      chr1     11-20       - ---      chr1     21-30       + |      <NA>
#>   [2]      chr1     21-30       + ---      chr2     51-60       + |         1
#>            id.2
#>       <integer>
#>   [1]         1
#>   [2]         2
#>   -------
#>   regions: 3 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

####################################################################
# 2. Match loops with tiled genomic bins
####################################################################

data(GM12878_HiCCUPS)
loops <- GM12878_HiCCUPS |> 
    pin_by('first') |> 
    anchor_center() |> 
    mutate(width1 = 500) |> 
    pin_by('second') |> 
    anchor_center() |> 
    mutate(width2 = 500)

genomic_bins <- GenomeInfoDb::getChromInfoFromUCSC(
    'hg19', assembled.molecules.only = TRUE, as.Seqinfo = TRUE
) |> 
    GenomicRanges::tileGenome(tilewidth = 10000) |> 
    unlist() |> 
    mutate(binID = seq_len(plyranges::n()))

annotate(loops, genomic_bins, by = 'binID') |> 
    select(starts_with('binID'))
#> GInteractions object with 10898 interactions and 2 metadata columns:
#>           seqnames1             ranges1 strand1     seqnames2
#>               <Rle>           <IRanges>   <Rle>         <Rle>
#>       [1]     chr10 100227250-100227749       * ---     chr10
#>       [2]     chr10 100227250-100227749       * ---     chr10
#>       [3]     chr10 101192250-101192749       * ---     chr10
#>       [4]     chr10 101194750-101195249       * ---     chr10
#>       [5]     chr10 101602250-101602749       * ---     chr10
#>       ...       ...                 ...     ... ...       ...
#>   [10894]      chrX     9847250-9847749       * ---      chrX
#>   [10895]      chrX     9847250-9847749       * ---      chrX
#>   [10896]      chrX     9962250-9962749       * ---      chrX
#>   [10897]      chrX   99764750-99765249       * ---      chrX
#>   [10898]      chrX   99942250-99942749       * ---      chrX
#>                       ranges2 strand2 |   binID.1   binID.2
#>                     <IRanges>   <Rle> | <integer> <integer>
#>       [1] 100422250-100422749       * |    178070    178090
#>       [2] 101007250-101007749       * |    178070    178148
#>       [3] 101372250-101372749       * |    178167    178185
#>       [4] 101474750-101475249       * |    178167    178195
#>       [5] 101807250-101807749       * |    178208    178228
#>       ...                 ...     ... .       ...       ...
#>   [10894]   10087250-10087749       * |    289111    289135
#>   [10895]     9962250-9962749       * |    289111    289123
#>   [10896]   10087250-10087749       * |    289123    289135
#>   [10897] 100024750-100025249       * |    298103    298129
#>   [10898] 100022250-100022749       * |    298121    298129
#>   -------
#>   regions: 16123 ranges and 0 metadata columns
#>   seqinfo: 23 sequences from an unspecified genome; no seqlengths

####################################################################
# 3. Annotate interactions by a set of regulatory elements
####################################################################

data(ce10_ARCC)
data(ce10_REs)
annotate(ce10_ARCC, ce10_REs, by = 'annot') |> 
   count(annot.1, annot.2) |> 
   as.data.frame() |> 
   arrange(desc(n))
#>                annot.1             annot.2    n
#> 1      coding_promoter     coding_promoter 5563
#> 2      coding_promoter   putative_enhancer 2546
#> 3    putative_enhancer     coding_promoter 2540
#> 4    putative_enhancer   putative_enhancer 1814
#> 5      coding_promoter unassigned_promoter  809
#> 6  unassigned_promoter     coding_promoter  753
#> 7    putative_enhancer unassigned_promoter  566
#> 8  unassigned_promoter   putative_enhancer  533
#> 9       non-coding_RNA     coding_promoter  247
#> 10     coding_promoter      non-coding_RNA  242
#> 11 unassigned_promoter unassigned_promoter  209
#> 12     coding_promoter       other_element  188
#> 13       other_element     coding_promoter  173
#> 14       other_element   putative_enhancer  115
#> 15      non-coding_RNA   putative_enhancer  110
#> 16   putative_enhancer       other_element  100
#> 17   putative_enhancer      non-coding_RNA   93
#> 18                <NA>     coding_promoter   54
#> 19 pseudogene_promoter     coding_promoter   45
#> 20   putative_enhancer                <NA>   44
#> 21      non-coding_RNA      non-coding_RNA   43
#> 22      non-coding_RNA unassigned_promoter   41
#> 23     coding_promoter pseudogene_promoter   40
#> 24     coding_promoter                <NA>   39
#> 25       other_element unassigned_promoter   36
#> 26       other_element       other_element   32
#> 27                <NA>   putative_enhancer   32
#> 28 unassigned_promoter      non-coding_RNA   25
#> 29 unassigned_promoter       other_element   25
#> 30       other_element      non-coding_RNA   19
#> 31                <NA>                <NA>   19
#> 32      non-coding_RNA       other_element   16
#> 33   putative_enhancer pseudogene_promoter   16
#> 34 unassigned_promoter                <NA>   13
#> 35 pseudogene_promoter   putative_enhancer   11
#> 36                <NA> unassigned_promoter   10
#> 37 pseudogene_promoter unassigned_promoter    9
#> 38       other_element pseudogene_promoter    6
#> 39 pseudogene_promoter pseudogene_promoter    5
#> 40       other_element                <NA>    3
#> 41 unassigned_promoter pseudogene_promoter    3
#> 42      non-coding_RNA                <NA>    2
#> 43 pseudogene_promoter      non-coding_RNA    2
#> 44 pseudogene_promoter       other_element    1
#> 45                <NA>      non-coding_RNA    1
#> 46                <NA>       other_element    1
```
