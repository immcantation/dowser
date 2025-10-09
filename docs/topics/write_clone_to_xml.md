**write_clone_to_xml** - *Takes an airr clone object and template and writes a BEAST2 XML file*

Description
--------------------

Takes an airr clone object and template and writes a BEAST2 XML file


Usage
--------------------
```
write_clone_to_xml(
clone,
file,
id,
time = NULL,
trait = NULL,
trait_data_type = NULL,
template = NULL,
mcmc_length = 1e+06,
log_every = 1000,
replacements = NULL,
include_germline_as_root = FALSE,
include_germline_as_tip = FALSE,
germline_range = c(-10000, 10000),
tree = NULL,
trait_list = NULL,
log_every_trait = 10,
tree_states = FALSE,
start_edge_length = 100,
start_date = NULL,
max_start_date = NULL,
...
)
```

Arguments
-------------------

clone
:   an `airrClone` object

file
:   output file path

id
:   unique identifer for this analysis

time
:   name of column representing sample time

trait
:   name of column representing a trait

trait_data_type
:   optional data type for the trait

template
:   XML template

mcmc_length
:   number of MCMC iterations

log_every
:   frequency of states logged. `auto` will divide mcmc_length by log_target

replacements
:   list of additional replacements to make in the template

include_germline_as_root
:   include germline in analysis as root?

include_germline_as_tip
:   include germline in analysis as tip?

germline_range
:   possible date range of germline

tree
:   starting tree, either a phylo object or a newick string

trait_list
:   list of all possible trait values

log_every_trait
:   frequency of trait states logged relative to log_every

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

File path of the written XML file









