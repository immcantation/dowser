**readBEAST** - *Reads in a BEAST output directory*

Description
--------------------

`readBEAST` Reads in data from BEAST output directory


Usage
--------------------
```
readBEAST(
clones,
dir,
id,
beast,
burnin = 10,
trait = NULL,
nproc = 1,
quiet = 0,
posterior = c("none", "all", "parameters", "trees_with_traits", "trees"),
asr = FALSE,
low_ram = TRUE,
trim_ids = FALSE
)
```

Arguments
-------------------

clones
:   either a tibble (getTrees) or list of `airrClone` objects

dir
:   directory where BEAST output files have been placed.

id
:   unique identifer for this analysis

beast
:   location of beast binary directory (beast/bin)

burnin
:   percent of initial tree samples to discard (default 10)

trait
:   Trait column used

nproc
:   Number of cores for parallelization. Uses at most 1 core per tree.

quiet
:   amount of rubbish to print to console

posterior
:   Read un full distribution of parameters and trees? Can be "none" to just have
summary objects, "all" to have parameters, trees, and trees_with_traits, or a vector with the desired
combination of "parameters", "trees_with_traits", and "trees".

asr
:   Log ancestral sequences?

low_ram
:   run with less memory (slightly slower)

trim_ids
:   remove last _ group from tips?




Value
-------------------

If data is a tibble, then the input clones tibble with additional columns for 
trees and parameter estimates given the specified burnin. If input is just a 
list of airrClone objects, it will return the corresponding list of trees
given the burnin









