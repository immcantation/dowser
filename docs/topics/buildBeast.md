**buildBeast** - *Read in a directory from a BEAST run. Runs treeannotator and loganalyser.*

Description
--------------------

Read in a directory from a BEAST run. Runs treeannotator and loganalyser.


Usage
--------------------
```
buildBeast(
data,
beast,
time,
template,
dir,
id,
mcmc_length = 1e+06,
resume_clones = NULL,
trait = NULL,
asr = FALSE,
full_posterior = FALSE,
log_every = "auto",
include_germline = TRUE,
nproc = 1,
quiet = 0,
burnin = 10,
low_ram = TRUE,
germline_range = c(-10000, 10000),
java = TRUE,
seed = NULL,
log_target = 10000,
trees = NULL,
tree_states = FALSE,
start_edge_length = 100,
start_date = NULL,
max_start_date = NULL,
...
)
```

Arguments
-------------------

data
:   a list of `airrClone` objects

beast
:   location of beast binary directory (beast/bin)

time
:   Name of sample time column

template
:   XML template

dir
:   directory where temporary files will be placed.

id
:   unique identifer for this analysis

mcmc_length
:   Number of MCMC iterations

resume_clones
:   Clones to resume for mcmc_length more iterations

trait
:   Trait coolumn used

asr
:   Log ancestral sequences?

full_posterior
:   Read un full distribution of parameters and trees?

log_every
:   Frequency of states logged. `auto` will divide mcmc_length by log_target

include_germline
:   Include germline in analysis?

nproc
:   Number of cores for parallelization. Uses 1 core per tree.

quiet
:   amount of rubbish to print to console

burnin
:   Burnin percent (default 10)

low_ram
:   run with less memory (slower)

germline_range
:   Possible date range of germline tip

java
:   Use the -java flag for BEAST run

seed
:   Used for the -seed option for BEASTrun

log_target
:   Target number of samples over mcmc_length

trees
:   optional list of starting trees, either phylo objects or newick strings

tree_states
:   Use `states` vector for starting tree

start_edge_length
:   edge length to use for all branches in starting tree

start_date
:   Starting date of time tree if desired

max_start_date
:   Maximum starting date of time tree if desired

...
:   Additional arguments for XML writing functions




Value
-------------------

The input clones tibble with an additional column for the bootstrap replicate trees.




See also
-------------------

[getTimeTrees](getTimeTrees.md)






