**readFromLastCheckpoint** - *Read the results of a run as they were at its last checkpoint*

Description
--------------------

Use if a `getTimeTreesIterate` run stopped and you want
the results as they were at the end of each clone's last finished round.


Usage
--------------------
```
readFromLastCheckpoint(
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
trim_ids = FALSE,
ess_cutoff = 200,
ignore = c("traitfrequencies")
)
```

Arguments
-------------------

clones
:   tibble of `airrClone` objects, as in `getTimeTreesIterate`

dir
:   directory of the analysis

id
:   analysis id

beast
:   location of the beast binary directory (beast/bin)

burnin
:   percent of initial samples to discard (default 10)

trait
:   trait column that was used, if any

nproc
:   number of cores for parallelization. Uses at most 1 core per clone.

quiet
:   amount of rubbish to print to console

posterior
:   what to read of the full posterior: "none" (default), "all", or a vector with "parameters", "trees_with_traits" and "trees". See `readBEAST`.

asr
:   were ancestral sequences logged?

low_ram
:   run with less memory (slightly slower)

trim_ids
:   remove last _ group from tips?

ess_cutoff
:   minimum ESS, for the column below_ESS

ignore
:   vector of parameters to ignore for the column below_ESS




Value
-------------------

Tibble of clones that have a good checkpoint, with additional columns
`trees`, `parameters`, `below_ESS` and `checkpoint_round` (number of rounds finished for that clone). A warning lists the clones that have no good checkpoint.




See also
-------------------

[getTimeTreesIterate](getTimeTreesIterate.md), [readBEAST](readBEAST.md)






