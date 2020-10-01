# Build trees

Dowser offer multiple ways to build B cell phylogenetic trees. These differ by the method used to estimate tree topology and branch lengths (e.g. maximum parsimony and maximum likelihood) and implementation (PHYLIP or R packages ape and phangorn).

Before trees can be built, B cell sequences must be separated into clonal clusters, and had their clonal germline sequences reconstructed. Default settings assume input data is in AIRR TSV format, though column names may be specified using function arguments.

## Format clones

Before trees can be built, data must be formatted into a data table of AIRR clone objects. This is accomplished using the formatClones function. This function will:

1. Change non-nucleotide characters to `N` characters.
2. By default, collapse sequences that are either identical or differ only by ambiguous characters. 
3. Sequences will not be collapsed if they differ by columns specified in the "traits" option, or if the `collapse` option is set to `FALSE`.
4. Include data columns specified by `num_fields` or `text_fields`.
5. Remove uninformative sequence sites in which all sequences have `N` characters.

The output of this function is a tibble in which each row is a clone, ordered by the number of sequences. The column `data` contains airrClone objects with the clonal sequence alignments. Other columns contain information about the clone, and can be specified using the `columns` argument.

```r
# load example AIRR tsv data
data(ExampleDb)

# Process example data using default settings
clones = formatClones(ExampleDb)

clones
## A tibble: 99 x 4
#   clone_id data       locus  seqs
#      <dbl> <list>     <chr> <int>
# 1     3128 <airrClon> N        43
# 2     3100 <airrClon> N        32
# 3     3141 <airrClon> N        17
# 4     3170 <airrClon> N        16


# Process example data keeping samples from different times
# distinct, adding duplicate_count among collapsed sequences,
# and show the d_call for each clone in the tibble.
clones = formatClones(ExampleDb, traits=c("sample_id"),
    num_fields=c("duplicate_count"), columns=c("d_call"))

clones
## A tibble: 100 x 5
#   clone_id data       locus  seqs d_call              
#      <dbl> <list>     <chr> <int> <chr>               
# 1     3128 <airrClon> N        43 Homsap IGHD6-19*01 F
# 2     3100 <airrClon> N        32 Homsap IGHD6-13*01 F
# 3     3141 <airrClon> N        17 Homsap IGHD6-19*01 F
# 4     3170 <airrClon> N        16 Homsap IGHD6-19*01 F
```

## Build trees using maximum parsimony

A common way to build B cell lineage trees is the find the tree topology that minimizes the number of mutations needed along the tree (i.e. is the most parsimonious). Branch lengths can then be estimated as the number of mutations per site between each node in the tree.

Maximum parsimony trees can be built with the `getTrees` function, which by default uses the `pratchet` maximum parsimony function in the `phangorn` phylogenetics package.

Maximum parsimony trees can also be built using the PHYLIP function `dnapars`. To do this, the `build` option needs to be set as `dnapars` and the path to the `dnapars` executable needs to be specified in the `exec` option.

In all cases, the output is the same tibble as the input, but with a `trees` column containing an R `ape::phylo` object for each clone.

```r

# Build trees using the pratchet maximum parsimony function in phangnorn
clones = getTrees(clones)

clones
## A tibble: 100 x 6
#   clone_id data       locus  seqs d_call               trees  
#      <dbl> <list>     <chr> <int> <chr>                <list> 
# 1     3128 <airrClon> N        43 Homsap IGHD6-19*01 F <phylo>
# 2     3100 <airrClon> N        32 Homsap IGHD6-13*01 F <phylo>
# 3     3141 <airrClon> N        17 Homsap IGHD6-19*01 F <phylo>
# 4     3170 <airrClon> N        16 Homsap IGHD6-19*01 F <phylo>

# Build trees using dnapars.
# exec here is set to dnapars position in the Docker image.
clones = getTrees(clones, build="dnapars", exec="/usr/local/bin/dnapars")

clones
## A tibble: 100 x 6
#   clone_id data       locus  seqs d_call               trees  
#      <dbl> <list>     <chr> <int> <chr>                <list> 
# 1     3128 <airrClon> N        43 Homsap IGHD6-19*01 F <phylo>
# 2     3100 <airrClon> N        32 Homsap IGHD6-13*01 F <phylo>
# 3     3141 <airrClon> N        17 Homsap IGHD6-19*01 F <phylo>
# 4     3170 <airrClon> N        16 Homsap IGHD6-19*01 F <phylo>

```

## Build trees using maximum likelihood

Coming soon!

