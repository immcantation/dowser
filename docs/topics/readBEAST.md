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
full_posterior = FALSE,
asr = FALSE,
low_ram = TRUE
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
:   Trait coolumn used

nproc
:   Number of cores for parallelization. Uses at most 1 core per tree.

quiet
:   amount of rubbish to print to console

full_posterior
:   Read un full distribution of parameters and trees?

asr
:   Log ancestral sequences?

low_ram
:   run with less memory (slightly slower)




Value
-------------------

If data is a tibble, then the input clones tibble with additional columns for 
trees and parameter estimates given the specified burnin. If input is just a 
list of airrClone objects, it will return the corresponding list of trees
given the burnin









