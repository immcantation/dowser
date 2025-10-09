**write_clones_to_xmls** - *Wrapper to write multiple clones to XML files*

Description
--------------------

Wrapper to write multiple clones to XML files


Usage
--------------------
```
write_clones_to_xmls(
data,
id,
trees = NULL,
time = NULL,
trait = NULL,
template = NULL,
outfile = NULL,
replacements = NULL,
trait_list = NULL,
mcmc_length = 1e+06,
log_every = 1000,
include_germline_as_root = FALSE,
include_germline_as_tip = FALSE,
germline_range = c(-10000, 10000),
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

id
:   identifer for this analysis

trees
:   optional list of starting trees, either phylo objects or newick strings

time
:   name of column representing sample time

trait
:   name of column representing a trait

template
:   XML template

outfile
:   output file path prefix

replacements
:   list of additional replacements to make in the template

trait_list
:   list of all possible trait values

mcmc_length
:   number of MCMC iterations

log_every
:   frequency of states logged. `auto` will divide mcmc_length by log_target

include_germline_as_root
:   include germline in analysis as root?

include_germline_as_tip
:   include germline in analysis as tip?

germline_range
:   possible date range of germline

tree_states
:   use states in the starting tree?

start_edge_length
:   edge length to use for all branches in starting tree

start_date
:   starting date to use as prior, in forward time

max_start_date
:   max starting date to use as prior, in forward time

...
:   additional arguments for XML writing functions




Value
-------------------

File paths of the written XML files









