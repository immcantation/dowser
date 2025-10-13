# Using data from non-B cells

While originally designed for B cells, Dowser also supports phylogenetic inference for non B cells, especially cells evolving from a known ancestral sequence, such as tumor lineages.

If sequences are from a single lineage, the only requirement for non-B cell data is that the sequences supplied to `formatClones` are aligned and in a data.frame with a column for sequences and a column for sequence IDs. If from multiple lineages, they can be deliminated using the `clone_id` column.

In the code block below, we show how trees can be built using data with only sequence IDs, sequences, and germline sequences.


``` r
library(dowser)

data(ExampleAirr)

ExampleAirr <- dplyr::select(dplyr::filter(ExampleAirr, clone_id=="3128"),
 sequence_id, sequence_alignment, germline_alignment)

clones <- formatClones(ExampleAirr, germ="germline_alignment")
```

```
## [1] "v_call, j_call, and junc_len not found in data. Using non B cell mode\n. Setting use_regions to FALSE."
```

``` r
trees <- getTrees(clones)
```

Note that if specified `v_call`, `j_call`, and `junction_length` columns are not found in the input data.frame, the options `use_regions` will be set to false, as it is only for BCR sequences. If not already present, the `clone_id` and `locus` columns will be added to the dataframe with values 0 and "N", respectively.

When using `getTimeTrees` or `getTimeTreesIterate`, a meaninful germline is not required. Instead, you can set a `germline_alignment` which is series of N nucleotides the same length as the sequences in the sequence_alignment column, and set `include_germline=FALSE.`







