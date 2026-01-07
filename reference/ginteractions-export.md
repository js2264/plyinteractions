# Export GInteractions as `bedpe` or `pairs` files

`write_*` functions are provided to export a GInteractions object into
these two file formats. See 4DN documentation
(https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md)
and UCSC documentation
(https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format)
for more details.

## Usage

``` r
write_bedpe(x, file, scores = NULL)

write_pairs(x, file, seqlengths = Seqinfo::seqlengths(x))
```

## Arguments

- x:

  a GInteractions object.

- file:

  path to a `.bedpe` or `.pairs` file to save the genomic interactions.

- scores:

  Name of column to extract scores from.

- seqlengths:

  Named vector indicating the chromosome sizes.

## Value

TRUE

## Examples

``` r
gi <- read.table(text = "
chr1 100 200 chr1 5000 5100 bedpe_example1 30 + -
chr1 1000 5000 chr1 3000 3800 bedpe_example2 100 + -",
col.names = c(
  "seqnames1", "start1", "end1", 
  "seqnames2", "start2", "end2", "name", "score", "strand1", "strand2")
) |> as_ginteractions()

write_bedpe(gi, 'gi.bedpe')
#> [1] TRUE
write_pairs(gi, 'gi.pairs')
#> No `seqlengths` provided. The `chromsizes` will be inferred from the interactions and will most likely be inaccurate.
#> [1] TRUE
```
