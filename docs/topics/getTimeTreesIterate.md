**getTimeTreesIterate** - *Iteratively resume getTimeTrees until convergence, as defined by 
all parameters (except those in `ignore` vector) having ESS 
greater than or equal to the specified ess_cutoff*

Description
--------------------

`getTimeTreesIterate` Iteratively resume getTimeTrees til convergence.


Usage
--------------------
```
getTimeTreesIterate(
clones,
template,
beast,
dir,
id,
time,
iterations = 10,
ess_cutoff = 200,
burnin = 10,
ignore = c("traitfrequencies"),
continue = FALSE,
posterior = c("none", "all", "parameters", "trees_with_traits", "trees"),
trait = NULL,
nproc = 1,
quiet = 0,
asr = FALSE,
low_ram = TRUE,
trim_ids = FALSE,
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

iterations
:   Maximum number of times to resume MCMC chain

ess_cutoff
:   Minimum number of ESS for all parameters

burnin
:   Burnin percent (default 10)

ignore
:   Vector of parameters to ignore for ESS calculation

continue
:   If TRUE, will check for iteration folder and resume from last iteration if found (default FALSE)

posterior
:   Read un full distribution of parameters and trees? Can be "none" to just have
summary objects, "all" to have parameters, trees, and trees_with_traits, or a vector with the desired
combination of "parameters", "trees_with_traits", and "trees".

trait
:   Trait column to be used

nproc
:   Number of cores for parallelization. At most 1 core/tree can be used.

quiet
:   quiet notifications if > 0

asr
:   Log ancestral sequences?

low_ram
:   run with less memory (slightly slower)

trim_ids
:   remove last _ group from tips?

...
:   Additional arguments for getTimeTrees




Value
-------------------

A tibble of `tidytree` and `airrClone` objects.


Details
-------------------

For a full list of options (of which there are many more), see `getTimeTrees`

For examples and vignettes, see https://dowser.readthedocs.io









