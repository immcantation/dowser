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
posterior = c("none", "all", "parameters", "trees_with_traits", "trees"),
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
germline_trait_value = "?",
trait_list = NULL,
rootfreqs_class = NULL,
xml_check = NULL,
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
:   Number of MCMC steps

resume_clones
:   Clones to resume for `mcmc_length` more steps

trait
:   Trait column used

asr
:   Log ancestral sequences?

posterior
:   Read un full distribution of parameters and trees? Can be "none" to just have
summary objects, "all" to have parameters, trees, and trees_with_traits, or a vector with the desired
combination of "parameters", "trees_with_traits", and "trees".

log_every
:   Frequency of states logged. `auto` will divide mcmc_length by log_target

include_germline
:   Include germline in analysis?

nproc
:   Number of cores for parallelization. Uses at most 1 core per tree.

quiet
:   Amount of rubbish to print to console

burnin
:   Burnin percent (default 10)

low_ram
:   run with less memory (slightly slower)

germline_range
:   Possible date range of germline tip

java
:   Use the -java flag for BEAST run

seed
:   Use specified seeed for the -seed option for BEAST

log_target
:   Target number of samples over `mcmc_length`

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

germline_trait_value
:   trait value for germline, default '?' for ambiguous

trait_list
:   optional vector of all types of the trait. If NULL, it is made from the values in `data`.

rootfreqs_class
:   BEAST class for the sequence likelihood when the germline is the root
(`include_germline=TRUE` and the template has `${ROOTFREQS}`). `NULL` (default)
runs `beast -version` once and uses
`tyche.evolution.likelihood.RootFreqsTreeLikelihood` for TyCHE 1.0.1 or later,
otherwise `rootfreqs.TreeLikelihood` (needs the `rootfreqs` BEAST package).
If given, BEAST is not run to find the version.

xml_check
:   optional function called with the paths of the XML files after they are
written and before BEAST starts. If it stops, BEAST is not started. 
getTimeTreesIterate uses this to get the XML of a clone without running BEAST.

...
:   Additional arguments for XML writing functions




Value
-------------------

The input clones tibble with an additional column for the bootstrap replicate trees.




See also
-------------------

[getTimeTrees](getTimeTrees.md)






