**findSwitches** - *Create a bootstrap distribution for clone sequence alignments, and estimate 
trees for each bootstrap replicate.*

Description
--------------------

`findSwitches` Phylogenetic bootstrap function.


Usage
--------------------
```
findSwitches(
clones,
permutations,
trait,
igphyml,
fixtrees = FALSE,
downsample = TRUE,
tip_switch = 20,
nproc = 1,
dir = NULL,
id = NULL,
modelfile = NULL,
build = "pratchet",
exec = NULL,
quiet = 0,
rm_temp = TRUE,
palette = NULL,
resolve = 2,
rep = NULL,
keeptrees = FALSE,
lfile = NULL,
seq = NULL,
boot_part = "locus",
force_resolve = FALSE,
...
)
```

Arguments
-------------------

clones
:   tibble `airrClone` objects, the output of 
[formatClones](formatClones.md)

permutations
:   number of bootstrap replicates to perform

trait
:   trait to use for parsimony models

igphyml
:   location of igphyml executible

fixtrees
:   keep tree topologies fixed?
(bootstrapping will not be perfomed)

downsample
:   downsample clones to have a maximum specified tip/switch ratio?

tip_switch
:   maximum allowed tip/switch ratio if downsample=TRUE

nproc
:   number of cores to parallelize computations

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

quiet
:   amount of rubbish to print to console

rm_temp
:   remove temporary files (default=TRUE)

palette
:   a named vector specifying colors for each state

resolve
:   how should polytomies be resolved? 
0=none, 1=max parsminy, 2=max ambiguity + polytomy skipping,
3=max ambiguity

rep
:   current bootstrap replicate (experimental)

keeptrees
:   keep trees estimated from bootstrap replicates? (TRUE)

lfile
:   lineage file input to igphyml if desired (experimental)

seq
:   column name containing sequence information

boot_part
:   is  "locus" bootstrap columns for each locus separately

force_resolve
:   continue even if polytomy resolution fails?

...
:   additional arguments to be passed to tree building program




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
(as well as `trait`), the returned object 
will contain a `tibble` named "switches" containing switch count 
information. This object can be passed to [testSP](testSP.md) and other functions 
to perform parsimony based trait value tests. 

Trait values cannot contain values N, UCA, or NTIP. These are reserved for 
use by test statistic functions.



Examples
-------------------

```R
### Not run:
data(ExampleAirr)
# ExampleAirr$sample_id <- sample(ExampleAirr$sample_id)
# clones <- formatClones(ExampleAirr, trait="sample_id")
# 
# igphyml <- "~/apps/igphyml/src/igphyml"
# btrees <- findSwitches(clones[1:2,], permutations=10, nproc=1,
# igphyml=igphyml, trait="sample_id")
# plotTrees(btrees$trees[[4]])[[1]]
# testPS(btrees$switches)
```



See also
-------------------

Uses output from [formatClones](formatClones.md) with similar arguments to 
[getTrees](getTrees.md). Output can be visualized with [plotTrees](plotTrees.md), and tested
with [testPS](testPS.md), [testSC](testSC.md), and [testSP](testSP.md).






