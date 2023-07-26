# Building trees with paired heavy and light chain data

With the advances in sequencing, single cell datasets can now pair heavy and light chain sequences. To best incorporate the paired data into BCR phylogenetics, Dowser offers a way to identify clone_subgroups. This step is done before formatting clones and building trees and allows for accurate phylogenetic construction of BCRs with paired heavy and light chain data. 

## Resolve light chains 

To resolve the light chains within a clone, use the resolveLightChains function. This function will:

1. Pair heavy and light chains together by their cell_id and identify which `clone_subgroup` each pair belongs to. 
2. Assign heavy chains without an associated light chain to the subgroup containing the most similar paired heavy chain.

The output of this function is a tibble in which each row is a different sequence, with all of the previously included data along with a few more columns. The column `clone_subgroup` contains the subgroup assignment for that sequence within a given clone, with 1 being the largest. `clone_subgroup_id` combines the `clone_id` variable and the `clone_subgroup` variable by a "_". `vj_cell` combines the `vj_gene` and `vj_alt_cell` columns by a ",". 


```r
library(dowser)
library(ggtree)
# load example tsv data
data("ExampleMixedDb")

# find the clone subgroups 
ExampleMixedDb <- resolveLightChains(ExampleMixedDb)
print(ExampleMixedDb$clone_subgroup)
# [1] 1 1 2 1 1 1 1 1 1 1 1 1 1 2 1 1 2 1 1 1 1 1 1 1 1 1 1
```


```r
# run createGermlines -- this has already been done on this data
#ExampleMixedDb <- createGermlines(ExampleMixedDb, clone = clone_subgroup_id, nproc = 1)
```

## Format clones

The next step is to convert the data into airrClone objects that can be used for tree building. As with heavy chain sequences, this is done using the formatClones function. To build trees with paired heavy and light chains, specify `chain` = "HL". This will concatenate the paired heavy and light chains into a single sequence alignment. To only use the heavy chain, simply leave chain = "H", the default. For more information on formatClones see the [Building Trees Vignette](Building-Trees-Vignette.md).


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

For details on each of the different methods, including the specifics about different partition models, see the [Building Trees Vignette](Building-Trees-Vignette.md).


```r
# Building maximum likelihood trees with multiple partitions using IgPhyML.
# exec here is set to IgPhyML position in the Docker image.
clones <- getTrees(clones[1,], build="igphyml", nproc=1, partition="hl",
                   exec="/usr/local/share/igphyml/src/igphyml")
```



```r
plotTrees(clones)[[1]]+geom_tiplab()
```



```r
library(dowser)
library(ggtree)
# Load data instead of running phylip
data(hlClone)
plotTrees(hlClone)[[1]]+geom_tiplab()
```

![plot of chunk Resolve-Light-Chains-Vignette-6](figure/Resolve-Light-Chains-Vignette-6-1.png)

Building maximum likelihood trees with multiple partitions using *RAxML* instead, which is similar to partition = "hl" in IgPhyML. RAxML does not estimate omega values, so this parameter is not necessary. 

```r
# exec here is set to RAxML position in the Docker image.
clones = getTrees(clones, build="raxml", 
    exec="/usr/local/bin/raxml-ng", nproc=1, partition=TRUE)

print(clones)
## A tibble: 2 x 7
#  clone_id data       locus  seqs subject_id trees        parameters       
#     <dbl> <list>     <chr> <int> <chr>      <named list> <named list>     
#1     3170 <airrClon> N        13 Subject_1  <phylo>      <named list [13]>
#2     3184 <airrClon> N        12 Subject_1  <phylo>      <named list [13]>i
```


```r
plotTrees(clones)[[1]]+geom_tiplab()
```


```r
library(dowser)
library(ggtree)
# Load data instead of running phylip
data(hlClone)
# change the tree names
hlClone$ig_tree <- hlClone$trees
hlClone$trees <- hlClone$raxml_tree
plotTrees(hlClone)[[1]]+geom_tiplab()
```

![plot of chunk Resolve-Light-Chains-Vignette-9](figure/Resolve-Light-Chains-Vignette-9-1.png)
