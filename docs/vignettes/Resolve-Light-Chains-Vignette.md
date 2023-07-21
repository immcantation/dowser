# Resolve Light Chains

With the advances in sequencing, single cell datasets can now track which light chain associates with which heavy chain. To best incorporate the paired data into use in BCR phylogenetics, Dowser offers a way to identify clone_subgroups. This step is done before formatting clones and building trees and allows for accurate phylogenetic construction of BCRs with paired heavy and light chain data. 

## Resolve light chains 

To resolve the light chains within a clone, use the resolveLightChains function. This function will:

1. Pair heavy and light chains together by their cell_id and identify which `clone_subgroup` each pair belongs to. 
2. Assign heavy chains without an associating light chain a clone_subgroup to the most similar paired heavy chain within the same clone.

The output of this function is a tibble in which each row is a different sequence, with all of the previously included data along with a few more columns. The column `clone_subgroup` contains the subgroup assignment for that sequence within a given clone. `vj_clone` combines the `clone_id` variable and the `clone_subgroup` variable by a "_". `vj_cell` combines the `vj_gene` and `vj_alt_cell` columns by a ",". 


```r
library(dowser)
# load example tsv data
data("ExampleMixedDb")

# find the clone subgroups 
ExampleMixedDb <- resolveLightChains(ExampleMixedDb)
print(ExampleMixedDb$clone_subgroup)
# [1] 1 1 2 1 1 1 1 1 1 1 1 1 1 2 1 1 2 1 1 1 1 1 1 1 1 1 1
```


```r
# run createGermlines -- this has already been done on this data
#ExampleMixedDb <- createGermlines(ExampleMixedDb, nproc = 1)
```

## Format clones

With the `clone_subgroup` now calculated, the data can now be formatted into a data table of AIRR clone objects. This is done using formatClones and specifying the `chain` variable to be "HL". This concatenates the paired heavy and light chains into a single sequence alignment. For more information on formatClones see the [Building Trees Vignette](Building-Trees-Vignette.md).


```r
# for paired data, ensure that the chain is HL. collapse = F here to get multiple clones due to the sequences being very similar. 
clones <- formatClones(ExampleMixedDb, chain="HL", nproc=1, collapse = F)
print(clones)
## A tibble: 3 x 4
#  clone_id     data       locus       seqs    
#     <dbl>    <list>     <chr>       <int> 
#1     1     <airrClon>   IGH,IGK       11 
#2     2     <airrClon>   IGH,IGK        2 
#3     4     <airrClon>   IGH,IGK        2 
```
## Building trees 

Trees can now be built with any of the methods that dowser supports. The broad categories are maximum parsimony, maximum likelihood, and partitioned maximum likelihood. However, the partitioned maximum likelihood models especially utilize paired heavy and light chain data. The maximum likelihood partition models will create separate scalars for heavy and light chain branch estimation, allowing for the most accurate branch length reconstitution. Here, branch lengths represent the number of mutations per site along both heavy and light chains. This is especially apparent when mixing bulk and single cell data. 

For specifics on each of the different methods, including the specifics about the different partitions, see the [Building Trees Vignette](Building-Trees-Vignette.md).


```r
# Building maximum likelihood trees with multiple partitions using IgPhyML.
# exec here is set to IgPhyML position in the Docker image.
clones <- getTrees(clones[1,], build="igphyml", nproc=1, omega="e,e", rates="0,1", partition="hl",
                   exec="/usr/local/share/igphyml/src/igphyml")
plotTrees(clones)
```
![]("figure/Resolve-Light-Chains-Vignette-1-1.png")
