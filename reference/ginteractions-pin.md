# Pin GInteractions by anchors set (anchors1 or anchors2).

Pin GInteractions by anchors set (anchors1 or anchors2).

## Usage

``` r
pin(x, anchors)

pin_by(x, anchors)

pinned_anchors(x)

unpin(x)

# S4 method for class 'GroupedGInteractions,character'
pin(x, anchors)

# S4 method for class 'GroupedGInteractions,numeric'
pin(x, anchors)

# S4 method for class 'GInteractions,character'
pin(x, anchors)

# S4 method for class 'GInteractions,numeric'
pin(x, anchors)

# S4 method for class 'PinnedGInteractions,missing'
pin(x, anchors)

# S4 method for class 'PinnedGInteractions,character'
pin(x, anchors)

# S4 method for class 'PinnedGInteractions,numeric'
pin(x, anchors)

# S4 method for class 'AnchoredPinnedGInteractions,character'
pin(x, anchors)

# S4 method for class 'AnchoredPinnedGInteractions,numeric'
pin(x, anchors)

pin_first(x)

pin_second(x)

pin_anchors1(x)

pin_anchors2(x)

# S4 method for class 'AnchoredPinnedGInteractions'
unpin(x)

# S4 method for class 'PinnedGInteractions'
unpin(x)

# S4 method for class 'GInteractions'
unpin(x)

# S4 method for class 'PinnedGInteractions'
pinned_anchors(x)

# S4 method for class 'AnchoredPinnedGInteractions'
pinned_anchors(x)
```

## Arguments

- x:

  a GInteractions object

- anchors:

  Anchors to pin on ("first" or "second")

## Value

- `pin_*` functions return a PinnedGInteractions object.

- `pin` returns a numerical value indicating which set of anchors is
  pinned.

- `unpin` removes the pinning of a PinnedGInteractions object.

- `pinned_anchors` returns an (Anchored)GenomicRanges object
  corresponding to the pinned anchors of a PinnedGInteractions object.

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
# 1. Pin by first anchors
####################################################################

gi |> pin_by("first")
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors1
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     11-20       * ---      chr1     21-30       * |         cis
#>   [2]      chr1     11-20       * ---      chr1     51-55       * |         cis
#>   [3]      chr1     11-30       * ---      chr1     51-55       * |         cis
#>   [4]      chr1     11-30       * ---      chr2     51-60       * |       trans
#>           score
#>       <numeric>
#>   [1]  0.718270
#>   [2]  0.241314
#>   [3]  0.547043
#>   [4]  0.834802
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_first()
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors1
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     11-20       * ---      chr1     21-30       * |         cis
#>   [2]      chr1     11-20       * ---      chr1     51-55       * |         cis
#>   [3]      chr1     11-30       * ---      chr1     51-55       * |         cis
#>   [4]      chr1     11-30       * ---      chr2     51-60       * |       trans
#>           score
#>       <numeric>
#>   [1]  0.718270
#>   [2]  0.241314
#>   [3]  0.547043
#>   [4]  0.834802
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_anchors1()
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors1
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     11-20       * ---      chr1     21-30       * |         cis
#>   [2]      chr1     11-20       * ---      chr1     51-55       * |         cis
#>   [3]      chr1     11-30       * ---      chr1     51-55       * |         cis
#>   [4]      chr1     11-30       * ---      chr2     51-60       * |       trans
#>           score
#>       <numeric>
#>   [1]  0.718270
#>   [2]  0.241314
#>   [3]  0.547043
#>   [4]  0.834802
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

####################################################################
# 2. Pin by second anchors
####################################################################

gi |> pin_by("second")
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors2
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     11-20       * ---      chr1     21-30       * |         cis
#>   [2]      chr1     11-20       * ---      chr1     51-55       * |         cis
#>   [3]      chr1     11-30       * ---      chr1     51-55       * |         cis
#>   [4]      chr1     11-30       * ---      chr2     51-60       * |       trans
#>           score
#>       <numeric>
#>   [1]  0.718270
#>   [2]  0.241314
#>   [3]  0.547043
#>   [4]  0.834802
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_second()
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors2
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     11-20       * ---      chr1     21-30       * |         cis
#>   [2]      chr1     11-20       * ---      chr1     51-55       * |         cis
#>   [3]      chr1     11-30       * ---      chr1     51-55       * |         cis
#>   [4]      chr1     11-30       * ---      chr2     51-60       * |       trans
#>           score
#>       <numeric>
#>   [1]  0.718270
#>   [2]  0.241314
#>   [3]  0.547043
#>   [4]  0.834802
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

gi |> pin_anchors2()
#> PinnedGInteractions object with 4 interactions and 2 metadata columns:
#> Pinned on: anchors2
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     11-20       * ---      chr1     21-30       * |         cis
#>   [2]      chr1     11-20       * ---      chr1     51-55       * |         cis
#>   [3]      chr1     11-30       * ---      chr1     51-55       * |         cis
#>   [4]      chr1     11-30       * ---      chr2     51-60       * |       trans
#>           score
#>       <numeric>
#>   [1]  0.718270
#>   [2]  0.241314
#>   [3]  0.547043
#>   [4]  0.834802
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

####################################################################
# 3. Unpin
####################################################################

gi |> pin("second") |> unpin()
#> GInteractions object with 4 interactions and 2 metadata columns:
#>       seqnames1   ranges1 strand1     seqnames2   ranges2 strand2 |        type
#>           <Rle> <IRanges>   <Rle>         <Rle> <IRanges>   <Rle> | <character>
#>   [1]      chr1     11-20       * ---      chr1     21-30       * |         cis
#>   [2]      chr1     11-20       * ---      chr1     51-55       * |         cis
#>   [3]      chr1     11-30       * ---      chr1     51-55       * |         cis
#>   [4]      chr1     11-30       * ---      chr2     51-60       * |       trans
#>           score
#>       <numeric>
#>   [1]  0.718270
#>   [2]  0.241314
#>   [3]  0.547043
#>   [4]  0.834802
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
