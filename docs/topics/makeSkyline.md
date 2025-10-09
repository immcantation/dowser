**makeSkyline** - *get values for Bayesian Skyline plot*

Description
--------------------

`makeSkyline`


Usage
--------------------
```
makeSkyline(
logfile,
treesfile,
burnin,
bins = 100,
youngest = 0,
clone_id = NULL,
max_height = c("min", "median", "mean", "max")
)
```

Arguments
-------------------

logfile
:   Beast log file

treesfile
:   BEAST trees file

burnin
:   Burnin percentage (1-100)

bins
:   number of bins for plotting

youngest
:   timepoint of the most recently tip sampled (if 0, backward time used)

clone_id
:   name of the clone being analyzed (if desired)




Value
-------------------

Bayesian Skyline values for given clone









