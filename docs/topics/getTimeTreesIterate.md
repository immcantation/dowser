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
skip_errors = FALSE,
trait_list = NULL,
backup = FALSE,
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
:   Maximum number of times to run MCMC chain for each clone. The first run counts as 1.

ess_cutoff
:   Minimum number of ESS for all parameters

burnin
:   Burnin percent (default 10)

ignore
:   Vector of parameters to ignore for ESS calculation

continue
:   If TRUE, continue an earlier run of the same `id` in `dir`
from saved checkpoints. (default FALSE) Warning: if the settings or
the data are not the same as in the earlier run, files may be overwritten, as if continue=FALSE.
See details.

posterior
:   Read un full distribution of parameters and trees? Can be "none" to just have
summary objects, "all" to have parameters, trees, and trees_with_traits, or a vector with the desired
combination of "parameters", "trees_with_traits", and "trees".

trait
:   Trait column to be used. "?" is treated as 
ambiguous (any type) and is not counted as a type.

nproc
:   Number of processes for parallelization. At most 1 process/tree can be used. Use nproc=1 for Windows.

quiet
:   quiet notifications if > 0

asr
:   Log ancestral sequences?

low_ram
:   run with less memory (slightly slower)

trim_ids
:   remove last _ group from tips?

skip_errors
:   If TRUE, a clone that fails is dropped with a warning and the other clones
continue. If FALSE (default), any failing clone stops the analysis.

trait_list
:   Optional vector of all types of the trait. If NULL, it is made from the values
in all clones, so that every clone uses the same types.

backup
:   If TRUE, a backup folder with the existing files of this analysis (BEAST files and
checkpoint folder) is made before the run starts. The files are moved if `continue`
is FALSE and copied if it is TRUE. See details. (default FALSE)

...
:   Additional arguments for getTimeTrees. If an argument is a function, it is
called with the one-row tibble of each clone, and the value that it returns
is used for that clone. A list in `trees` is cut to the tree of each clone.




Value
-------------------

A tibble of `tidytree` and `airrClone` objects. Clones that failed with `skip_errors=TRUE` are not in the tibble.


Details
-------------------

For a full list of options (of which there are many more), see `getTimeTrees`

Each clone runs as its own BEAST job. A clone with ESS below `ess_cutoff` is resumed (its MCMC
chain continues) until it converges or has run `iterations` times. For resumability in case of interruption, some intermediate and checkpoint files are saved. 

**continue=TRUE** is for resuming a job that was cancelled, or for running more `iterations`.
Repeat the original call, with the same `clones`, `dir`, `id` and settings, and add
`continue=TRUE` (with a larger `iterations` to run more). Clones with a checkpoint resume from
it, unless they already converged or reached `iterations`. Clones without one start from the
beginning.
You can also change `ess_cutoff` or `ignore`: convergence is checked again from the last
checkpoint of each clone, and clones that are now below the cutoff run again if they have run fewer than
`iterations` times.

Do not use `continue=TRUE` when changing settings (e.g. `seed`, `mcmc_length` or
`template`) or the data of an analysis that you want to keep. Each checkpoint has the XML that was
used for the clone, and the seed, the java flag and the BEAST path. If any of them is not the same in this
call, the checkpoint and the BEAST files of that clone are deleted, and the clone starts again from round 1.
The other clones continue, and a message lists the clones. Use a new `id` or `dir` instead. Do not
use if you have edited the BEAST files in `dir`.
Checkpoints from older versions of dowser cannot be used.

**backup=TRUE** makes the folder `backup_<id>_<date>_<time>` in `dir` before the run starts, with the
BEAST and checkpoint files of `id`. Use it if your job may be
restarted automatically (for example by a job scheduler), e.g. the restarted job starts a new run, you cancel it and
use `[restoreIterateBackup](restoreIterateBackup.md)` to put the files back, and then call the function again with
`continue=TRUE`. Use it also for insurance if `dir` may not be empty before this call.

If a clone fails and `skip_errors=FALSE`, running clones finish their current round, no new clone
starts, and then the function stops with the error. 

If you interrupt R (for example with Ctrl-C), BEAST processes that the workers started can keep running.
Check for them with the appropriate tools for your system.

For examples and vignettes, see https://dowser.readthedocs.io









