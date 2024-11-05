**bootstrapTrees** - *Deprecated! Please use findSwitches instead.*

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
fixtrees = FALSE,
quiet = 0,
rm_temp = TRUE,
palette = NULL,
resolve = 2,
rep = NULL,
keeptrees = TRUE,
lfile = NULL,
seq = NULL,
downsample = FALSE,
tip_switch = 20,
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
:   unique identifier for this analysis (required if 
`igphyml` or `dnapars` specified)

modelfile
:   file specifying parsimony model to use

build
:   program to use for tree building (phangorn, dnapars)

exec
:   location of desired phylogenetic executable

igphyml
:   location of igphyml executable if trait models desired

fixtrees
:   keep tree topologies fixed?
(bootstrapping will not be performed)

quiet
:   amount of rubbish to print to console

rm_temp
:   remove temporary files (default=TRUE)

palette
:   deprecated

resolve
:   how should polytomies be resolved? 
0=none, 1=max parsimony, 2=max ambiguity + polytomy skipping,
3=max ambiguity

rep
:   current bootstrap replicate (experimental)

keeptrees
:   keep trees estimated from bootstrap replicates? (TRUE)

lfile
:   lineage file input to igphyml if desired (experimental)

seq
:   column name containing sequence information

downsample
:   downsample clones to have a maximum specified tip/switch ratio?

tip_switch
:   maximum allowed tip/switch ratio if downsample=TRUE

boot_part
:   is  "locus" bootstrap columns for each locus separately

force_resolve
:   continue even if polytomy resolution fails?

...
:   additional arguments to be passed to tree building program




Value
-------------------

A list of trees and/or switch counts for each bootstrap replicate.









