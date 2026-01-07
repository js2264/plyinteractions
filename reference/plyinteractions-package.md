# plyinteractions: Extending tidyomics verbs to genomic interactions

plyinteractions verbs treat `GInteractions` objects as tabular data
using `dplyr`-like verbs. The functions and methods in `plyinteractions`
provide a grammatical approach to manipulate `GInteractions`, to
facilitate their integration in genomic analysis workflows.

plyinteractions is a dplyr-like API to the GInteractions infrastructure
in Bioconductor.

## Details

plyinteractions provides a consistent interface for importing and
wrangling genomic interactions from a variety of sources. The package
defines a grammar of genomic interactions manipulation through a set of
verbs. These verbs can be used to construct human-readable analysis
pipelines based on `GInteractions`.

- Group genomic interactions with [`group_by`](dplyr-group_by.md);

- Summarize grouped genomic interactions with
  [`summarize`](dplyr-summarize.md);

- Tally/count grouped genomic interactions with
  [`tally`](dplyr-count.md) and [`count`](dplyr-count.md);

- Modify genomic interactions with [`mutate`](dplyr-mutate.md);

- Subset genomic interactions with [`filter`](dplyr-filter.md) using
  [`<data-masking>`](https://rlang.r-lib.org/reference/args_data_masking.html)
  and logical expressions;

- Pick out any columns from the associated metadata with
  [`select`](dplyr-select.md) using [`<tidy-select>`
  arguments](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html);

- Subset using indices with [`slice`](dplyr-slice.md);

- Order genomic interactions with [`arrange`](dplyr-arrange.md) using
  categorical/numerical variables.

  For more details on the features of plyinteractions, read the
  vignette: `browseVignettes(package = "plyinteractions")`

## See also

Useful links:

- <https://github.com/js2264/plyinteractions>

- Report bugs at <https://github.com/js2264/plyinteractions/issues>

Useful links:

- <https://github.com/js2264/plyinteractions>

- Report bugs at <https://github.com/js2264/plyinteractions/issues>

## Author

Jacques Serizay

**Maintainer**: Jacques Serizay <jacquesserizay@gmail.com>
