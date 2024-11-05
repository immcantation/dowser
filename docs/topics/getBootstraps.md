**getBootstraps** - *Creates a bootstrap distribution for clone sequence alignments, and returns  
estimated trees for each bootstrap replicate as a nested list as a new input 
tibble column.*

Description
--------------------

`getBootstraps` Phylogenetic bootstrap function.


Usage
--------------------
```
getBootstraps(
clones,
bootstraps,
nproc = 1,
bootstrap_nodes = TRUE,
dir = NULL,
id = NULL,
build = "pratchet",
exec = NULL,
quiet = 0,
rm_temp = TRUE,
rep = NULL,
seq = NULL,
boot_part = "locus",
by_codon = TRUE,
starting_tree = FALSE,
switches = FALSE,
...
)
```

Arguments
-------------------

clones
:   tibble `airrClone` objects, the output of 
[formatClones](formatClones.md)

bootstraps
:   number of bootstrap replicates to perform

nproc
:   number of cores to parallelize computations

bootstrap_nodes
:   a logical if the the nodes for each tree in the trees
column (required) should report their bootstrap value

dir
:   directory where temporary files will be placed (required
if `igphyml` or `dnapars` specified)

id
:   unique identifier for this analysis (required if 
`igphyml` or `dnapars` specified)

build
:   program to use for tree building (phangorn, dnapars, igphyml)

exec
:   location of desired phylogenetic executable

quiet
:   amount of rubbish to print to console

rm_temp
:   remove temporary files (default=TRUE)

rep
:   current bootstrap replicate (experimental)

seq
:   column name containing sequence information

boot_part
:   is  "locus" bootstrap columns for each locus separately

by_codon
:   a logical if the user wants to bootstrap by codon or by 
nucleotide. Default (codon based bootstrapping) is TRUE.

starting_tree
:   An indicator to use the existing trees column as the starting trees for RAxML

switches
:   a logical indicator to allow findSwitches to do permutations.

...
:   additional arguments to be passed to tree building program




Value
-------------------

The input clones tibble with an additional column for the bootstrap replicate trees.









