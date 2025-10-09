**getSkylines** - *Make data frames for Bayesian skyline plots*

Description
--------------------

`makeSkylines`


Usage
--------------------
```
getSkylines(
clones,
dir,
id,
time,
burnin = 10,
bins = 100,
verbose = 0,
forward = TRUE,
nproc = 1,
max_height = c("min", "median", "mean", "max")
)
```

Arguments
-------------------

clones
:   clone tibble

dir
:   directory of BEAST trees file

time
:   name of time column

burnin
:   Burnin percent (default 10)

bins
:   number of bins for plotting

verbose
:   if 1, print name of clones

forward
:   plot in forward or (FALSE) backward time?

nproc
:   processors for parallelization (by clone)




Value
-------------------

Bayesian Skyline values for given clone


Details
-------------------

Burnin set from readBEAST or getTrees









