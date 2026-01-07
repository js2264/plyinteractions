# Manage GInteractions anchors with plyranges

Manage GInteractions anchors with plyranges

## Usage

``` r
# S3 method for class 'AnchoredPinnedGInteractions'
anchor(x)

# S3 method for class 'AnchoredPinnedGInteractions'
unanchor(x)

# S3 method for class 'PinnedGInteractions'
anchor_start(x)

# S3 method for class 'PinnedGInteractions'
anchor_end(x)

# S3 method for class 'PinnedGInteractions'
anchor_center(x)

# S3 method for class 'PinnedGInteractions'
anchor_3p(x)

# S3 method for class 'PinnedGInteractions'
anchor_5p(x)

# S3 method for class 'AnchoredPinnedGInteractions'
anchor_start(x)

# S3 method for class 'AnchoredPinnedGInteractions'
anchor_end(x)

# S3 method for class 'AnchoredPinnedGInteractions'
anchor_center(x)

# S3 method for class 'AnchoredPinnedGInteractions'
anchor_3p(x)

# S3 method for class 'AnchoredPinnedGInteractions'
anchor_5p(x)
```

## Arguments

- x:

  A PinnedGInteractions object

## Value

- `anchor_*` functions return an AnchoredPinnedGInteractions object.

- `anchor` returns a character string indicating where the pinned
  anchors are anchored at.

- `unanchor` removes the anchoring for a AnchoredPinnedGInteractions
  object.

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

gi
#> GInteractions object with 4 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1     11-20       + ---      chr1     21-30       + | 0.0280610
#>   [2]      chr1     11-20       + ---      chr1     51-55       + | 0.4659872
#>   [3]      chr1     11-30       - ---      chr1     51-55       - | 0.3900314
#>   [4]      chr1     11-30       - ---      chr2     51-60       - | 0.0200652
#>              type
#>       <character>
#>   [1]         cis
#>   [2]         cis
#>   [3]         cis
#>   [4]       trans
#>   -------
#>   regions: 6 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

####################################################################
# 1. Anchoring pinned genomic interactions with plyranges
####################################################################

gi |> pin_by("second") |> anchor_end()
#> AnchoredPinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors2 | Anchored by: end
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1     11-20       + ---      chr1     21-30       + | 0.0280610
#>   [2]      chr1     11-20       + ---      chr1     51-55       + | 0.4659872
#>   [3]      chr1     11-30       - ---      chr1     51-55       - | 0.3900314
#>   [4]      chr1     11-30       - ---      chr2     51-60       - | 0.0200652
#>              type
#>       <character>
#>   [1]         cis
#>   [2]         cis
#>   [3]         cis
#>   [4]       trans
#>   -------
#>   regions: 6 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_by("first") |> anchor_start()
#> AnchoredPinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors1 | Anchored by: start
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1     11-20       + ---      chr1     21-30       + | 0.0280610
#>   [2]      chr1     11-20       + ---      chr1     51-55       + | 0.4659872
#>   [3]      chr1     11-30       - ---      chr1     51-55       - | 0.3900314
#>   [4]      chr1     11-30       - ---      chr2     51-60       - | 0.0200652
#>              type
#>       <character>
#>   [1]         cis
#>   [2]         cis
#>   [3]         cis
#>   [4]       trans
#>   -------
#>   regions: 6 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> 
  pin_by("first") |> anchor_center() |> stretch(4) |> 
  pin_by("second") |> anchor_3p() |> stretch(-2)
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors2
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |     score
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <numeric>
#>   [1]      chr1      9-22       + ---      chr1     23-30       + | 0.0280610
#>   [2]      chr1      9-22       + ---      chr1     53-55       + | 0.4659872
#>   [3]      chr1      9-32       - ---      chr1     51-53       - | 0.3900314
#>   [4]      chr1      9-32       - ---      chr2     51-58       - | 0.0200652
#>              type
#>       <character>
#>   [1]         cis
#>   [2]         cis
#>   [3]         cis
#>   [4]       trans
#>   -------
#>   regions: 6 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
