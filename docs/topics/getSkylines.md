**getSkylines** - *Make data frames for Bayesian skyline plots*

Description
--------------------

`makeSkylines`


Usage
--------------------
```
getSkylines(
clones,
time,
bins = 100,
verbose = 0,
forward = TRUE,
nproc = 1,
max_height = c("min", "median", "mean", "max"),
exclude_germline = TRUE
)
```

Arguments
-------------------

clones
:   clone tibble

time
:   name of time column

bins
:   number of bins for plotting

verbose
:   if 1, print name of clones

forward
:   plot in forward or (FALSE) backward time?

nproc
:   processors for parallelization (by clone)

max_height
:   max height to use (min, median, mean, max)

exclude_germline
:   exclude germline from skyline plot? (For TyCHE GRTs)




Value
-------------------

Bayesian Skyline values for given clone


Details
-------------------

Clones must contain treedata objects with parameters_posterior and
trees_posterior. See `readBEAST` or `getTimeTrees` with posterior="all"









