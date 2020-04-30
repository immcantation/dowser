**bootstrapTrees** - *Create a bootstrap distribution for clone sequence alignments, and estimate 
trees for each bootstrap replicate.*

Description
--------------------

`bootstrapTrees` Phylogenetic bootstrap function.


Usage
--------------------
```
bootstrapTrees(
clones,
bootstraps,
nproc = 1,
trait = NULL,
dir = NULL,
id = NULL,
modelfile = NULL,
build = "pratchet",
exec = NULL,
igphyml = NULL,
trees = NULL,
quiet = 0,
rm_temp = TRUE,
palette = NULL,
resolve = 2,
rep = NULL,
keeptrees = TRUE,
lfile = NULL,
seq = "sequence"
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

trait
:   trait to use for parsimony models (required if 
`igphyml` specified)

dir
:   directory where temporary files will be placed (required
if `igphyml` or `dnapars` specified)

id
:   unique identifer for this analysis (required if 
`igphyml` or `dnapars` specified)

modelfile
:   file specifying parsimony model to use

build
:   program to use for tree building (phangorn, dnapars)

exec
:   location of desired phylogenetic executable

igphyml
:   location of igphyml executible if trait models desired

trees
:   tree topologies to use if aready available 
(bootstrapping will not be perfomed)

quiet
:   amount of rubbish to print to console

rm_temp
:   remove temporary files (default=TRUE)

palette
:   a named vector specifying colors for each state

resolve
:   how should polytomies be resolved?

rep
:   current bootstrap replicate (experimental)

keeptrees
:   keep trees estimated from bootstrap replicates? (TRUE)

lfile
:   lineage file input to igphyml if desired (experimental)

seq
:   column name containing sequence information




Value
-------------------

A list of trees and/or switch counts for each bootstrap replicate.


Details
-------------------

Tree building details are the same as [getTrees](getTrees.md). 
If `keeptrees=TRUE` (default) the returned object will contain a list 
named "trees" which contains a list of estimated tree objects for each 
bootstrap replicate. The object is structured like: 
trees[[<replicate>]][[<tree index>]]. If `igphyml` is specified 
(as well as `trait`, `temp`, and `id`), the returned object 
will contain a `tibble` named "switches" containing switch count 
information. This object can be passed to PStest and other functions 
to perform parsimony based trait value tests.



Examples
-------------------

```R
### Not run:
data(ExampleDb)
# ExampleDb$sample_id <- sample(ExampleDb$sample_id)
# clones <- formatClones(ExampleDb, trait="sample_id")
# 
# btrees <- bootstrapTrees(clones[1:2], bootstraps=100)
# plotTrees(btrees$trees[[4]][[1]])
# 
# igphyml <- "~/apps/igphyml/src/igphyml"
# btrees <- bootstrapTrees(clones[1:2], bootstraps=100, nproc=1,
# igphyml=igphyml, trait="sample_id", id="temp", dir="temp")
# plotTrees(btrees$trees[[4]][[1]])
# testPS(btrees$switches)
```



See also
-------------------

Uses output from [formatClones](formatClones.md) with similar arguments to 
[getTrees](getTrees.md). Output can be visualized with [plotTrees](plotTrees.md), and tested
with [testPS](testPS.md), [testSC](testSC.md), and [testSP](testSP.md).






