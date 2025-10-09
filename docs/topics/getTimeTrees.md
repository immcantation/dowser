**getTimeTrees** - *Estimate time trees by running BEAST on each clone
Applies XML `template` to each clone*

Description
--------------------

`getTimeTrees` Tree building function.


Usage
--------------------
```
getTimeTrees(
clones,
template,
beast,
dir,
id,
time,
mcmc_length = 3e+07,
log_every = "auto",
burnin = 10,
trait = NULL,
resume_clones = NULL,
nproc = 1,
quiet = 0,
rm_temp = FALSE,
include_germline = TRUE,
seq = "sequence",
germline_range = c(-10000, 10000),
java = TRUE,
seed = NULL,
log_target = 10000,
tree_states = FALSE,
trees = NULL,
...
)
```

Arguments
-------------------

clones
:   a tibble of `airrClone` objects, the output of
[formatClones](formatClones.md)

template
:   XML template

beast
:   location of beast binary directory (beast/bin)

dir
:   directory where temporary files will be placed.

id
:   unique identifer for this analysis

time
:   Name of sample time column

mcmc_length
:   Number of MCMC iterations

log_every
:   Frequency of states logged. "auto" will divide
mcmc_length by log_target

burnin
:   Burnin percent (default 10)

trait
:   Trait coolumn used

resume_clones
:   Clones to resume for mcmc_length more iterations

nproc
:   Number of cores for parallelization. Uses 1 core per tree.

quiet
:   amount of rubbish to print to console

rm_temp
:   remove temporary files (default=TRUE)

include_germline
:   Include germline in analysis?

seq
:   Sequence column in data

germline_range
:   Possible date range of germline tip

java
:   Use the -java flag for BEAST run

seed
:   Used for the -seed option for BEASTrun

log_target
:   Target number of samples over mcmc_length

tree_states
:   Use `states` vector for starting tree

trees
:   optional list of starting trees, either phylo objects or newick strings

...
:   Additional arguments passed to tree building programs




Value
-------------------

A list of `phylo` objects in the same order as `data`.


Details
-------------------

For examples and vignettes, see https://dowser.readthedocs.io




See also
-------------------

[getTrees](getTrees.md), [readBEAST](readBEAST.md)






