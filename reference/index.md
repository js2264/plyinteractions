# Package index

## About plyinteractions

- [`plyinteractions-package`](plyinteractions-package.md)
  [`plyinteractions`](plyinteractions-package.md) : plyinteractions:
  Extending tidyomics verbs to genomic interactions

## Constructors

- [`as_ginteractions()`](ginteractions-construct.md) : Construct a
  GInteractions object from a tibble, DataFrame or data.frame
- [`anchors1()`](ginteractions-getters.md)
  [`anchors2()`](ginteractions-getters.md)
  [`seqnames1()`](ginteractions-getters.md)
  [`seqnames2()`](ginteractions-getters.md)
  [`start1()`](ginteractions-getters.md)
  [`start2()`](ginteractions-getters.md)
  [`end1()`](ginteractions-getters.md)
  [`end2()`](ginteractions-getters.md)
  [`width1()`](ginteractions-getters.md)
  [`width2()`](ginteractions-getters.md)
  [`strand1()`](ginteractions-getters.md)
  [`strand2()`](ginteractions-getters.md)
  [`ranges1()`](ginteractions-getters.md)
  [`ranges2()`](ginteractions-getters.md)
  [`` `$`( ``*`<GInteractions>`*`)`](ginteractions-getters.md) :
  Enhanced GInteractions getters
- [`pair_granges()`](pair-granges.md) : Pairwise combination of a
  GRanges object

## `dplyr` core verbs

- [`arrange(`*`<GInteractions>`*`)`](dplyr-arrange.md) : Arrange a
  GInteractions by a column

- [`tally(`*`<GroupedGInteractions>`*`)`](dplyr-count.md)
  [`count(`*`<GroupedGInteractions>`*`)`](dplyr-count.md)
  [`count(`*`<GInteractions>`*`)`](dplyr-count.md) : Count or tally
  GInteractions per group

- [`filter(`*`<GInteractions>`*`)`](dplyr-filter.md) :

  Subset a GInteractions with tidyverse-like `filter`

- [`group_by(`*`<GInteractions>`*`)`](dplyr-group_by.md)
  [`group_by(`*`<DelegatingGInteractions>`*`)`](dplyr-group_by.md)
  [`ungroup(`*`<GroupedGInteractions>`*`)`](dplyr-group_by.md) : Group
  GInteractions by columns

- [`mutate(`*`<GInteractions>`*`)`](dplyr-mutate.md) : Mutate columns
  from a GInteractions object

- [`rename(`*`<GInteractions>`*`)`](dplyr-rename.md) :

  Rename columns from a GInteractions with tidyverse-like `rename`

- [`select(`*`<GInteractions>`*`)`](dplyr-select.md) : Select columns
  within GInteractions metadata columns

- [`slice(`*`<GInteractions>`*`)`](dplyr-slice.md) : Slice a
  GInteractions rows by their index

- [`summarise(`*`<GroupedGInteractions>`*`)`](dplyr-summarize.md)
  [`summarize(`*`<GroupedGInteractions>`*`)`](dplyr-summarize.md) :
  Summarize GInteractions per group

## `dplyr` group helpers

- [`group_data(`*`<GroupedGInteractions>`*`)`](group-group_data.md)
  [`group_keys(`*`<GroupedGInteractions>`*`)`](group-group_data.md)
  [`group_indices(`*`<GroupedGInteractions>`*`)`](group-group_data.md)
  [`group_vars(`*`<GInteractions>`*`)`](group-group_data.md)
  [`group_vars(`*`<GroupedGInteractions>`*`)`](group-group_data.md)
  [`groups(`*`<GroupedGInteractions>`*`)`](group-group_data.md)
  [`group_size(`*`<GroupedGInteractions>`*`)`](group-group_data.md)
  [`n_groups(`*`<GroupedGInteractions>`*`)`](group-group_data.md) :
  GInteractions grouping metadata

## `plyranges` verbs

- [`flank_downstream()`](plyinteractions-flank.md)
  [`flank_left()`](plyinteractions-flank.md)
  [`flank_upstream()`](plyinteractions-flank.md)
  [`flank_right()`](plyinteractions-flank.md) : Generate flanking
  regions from pinned anchors of a GInteractions object
- [`shift_downstream()`](plyinteractions-shift.md)
  [`shift_upstream()`](plyinteractions-shift.md)
  [`shift_right()`](plyinteractions-shift.md)
  [`shift_left()`](plyinteractions-shift.md) : Shift pinned anchors of a
  GInteractions object with plyinteractions
- [`stretch(`*`<AnchoredPinnedGInteractions>`*`)`](plyranges-stretch.md)
  [`stretch(`*`<PinnedGInteractions>`*`)`](plyranges-stretch.md) :
  Stretch pinned anchors of a GInteractions object with plyranges

## Overlapping GInteractions

- [`find_overlaps(`*`<PinnedGInteractions>`*`)`](ginteractions-find-overlaps.md)
  [`find_overlaps(`*`<GInteractions>`*`)`](ginteractions-find-overlaps.md)
  [`find_overlaps_directed(`*`<PinnedGInteractions>`*`)`](ginteractions-find-overlaps.md)
  [`find_overlaps_directed(`*`<GInteractions>`*`)`](ginteractions-find-overlaps.md)
  : Find overlaps between a query GInteractions and a GRanges
- [`count_overlaps(`*`<PinnedGInteractions>`*`)`](ginteractions-count-overlaps.md)
  [`count_overlaps(`*`<GInteractions>`*`)`](ginteractions-count-overlaps.md)
  [`count_overlaps_directed(`*`<PinnedGInteractions>`*`)`](ginteractions-count-overlaps.md)
  [`count_overlaps_directed(`*`<GInteractions>`*`)`](ginteractions-count-overlaps.md)
  : Count overlaps between a query GInteractions and a GRanges
- [`filter_by_overlaps(`*`<PinnedGInteractions>`*`)`](ginteractions-filter-overlaps.md)
  [`filter_by_overlaps(`*`<GInteractions>`*`)`](ginteractions-filter-overlaps.md)
  [`filter_by_non_overlaps(`*`<PinnedGInteractions>`*`)`](ginteractions-filter-overlaps.md)
  [`filter_by_non_overlaps(`*`<GInteractions>`*`)`](ginteractions-filter-overlaps.md)
  : Filter GInteractions overlapping with a GRanges
- [`join_overlap_left(`*`<PinnedGInteractions>`*`)`](ginteractions-join-overlap-left.md)
  [`join_overlap_left(`*`<GInteractions>`*`)`](ginteractions-join-overlap-left.md)
  [`join_overlap_left_directed(`*`<PinnedGInteractions>`*`)`](ginteractions-join-overlap-left.md)
  [`join_overlap_left_directed(`*`<GInteractions>`*`)`](ginteractions-join-overlap-left.md)
  : Join overlaps between a query GInteractions and a GRanges

## Enriching GInteractions

- [`annotate()`](ginteractions-annotate.md)
  [`annotate_directed()`](ginteractions-annotate.md) : Annotate both
  anchors of a GInteractions
- [`add_pairdist()`](add-pairdist.md) : Appends distance between
  interaction anchors

## Pinning GInteractions

- [`pin()`](ginteractions-pin.md) [`pin_by()`](ginteractions-pin.md)
  [`pinned_anchors()`](ginteractions-pin.md)
  [`unpin()`](ginteractions-pin.md)
  [`pin_first()`](ginteractions-pin.md)
  [`pin_second()`](ginteractions-pin.md)
  [`pin_anchors1()`](ginteractions-pin.md)
  [`pin_anchors2()`](ginteractions-pin.md) : Pin GInteractions by
  anchors set (anchors1 or anchors2).

## Anchoring GInteractions

- [`anchor(`*`<AnchoredPinnedGInteractions>`*`)`](ginteractions-anchor.md)
  [`unanchor(`*`<AnchoredPinnedGInteractions>`*`)`](ginteractions-anchor.md)
  [`anchor_start(`*`<PinnedGInteractions>`*`)`](ginteractions-anchor.md)
  [`anchor_end(`*`<PinnedGInteractions>`*`)`](ginteractions-anchor.md)
  [`anchor_center(`*`<PinnedGInteractions>`*`)`](ginteractions-anchor.md)
  [`anchor_3p(`*`<PinnedGInteractions>`*`)`](ginteractions-anchor.md)
  [`anchor_5p(`*`<PinnedGInteractions>`*`)`](ginteractions-anchor.md)
  [`anchor_start(`*`<AnchoredPinnedGInteractions>`*`)`](ginteractions-anchor.md)
  [`anchor_end(`*`<AnchoredPinnedGInteractions>`*`)`](ginteractions-anchor.md)
  [`anchor_center(`*`<AnchoredPinnedGInteractions>`*`)`](ginteractions-anchor.md)
  [`anchor_3p(`*`<AnchoredPinnedGInteractions>`*`)`](ginteractions-anchor.md)
  [`anchor_5p(`*`<AnchoredPinnedGInteractions>`*`)`](ginteractions-anchor.md)
  : Manage GInteractions anchors with plyranges
- [`replace_anchors()`](replace_anchors.md) : Replace anchors of a
  GInteractions

## Exporting GInteractions

- [`write_bedpe()`](ginteractions-export.md)
  [`write_pairs()`](ginteractions-export.md) :

  Export GInteractions as `bedpe` or `pairs` files

## Misc

- [`show-GInteractions`](show-GInteractions.md)
  [`show,GInteractions-method`](show-GInteractions.md)
  [`show,AnchoredPinnedGInteractions-method`](show-GInteractions.md)
  [`show,GroupedGInteractions-method`](show-GInteractions.md)
  [`show,PinnedGInteractions-method`](show-GInteractions.md) :

  show method for `GInteractions` objects

## Toy datasets

- [`GM12878_HiCCUPS`](plyinteractions-data.md) : Data files provided in
  the plyinteractions package
