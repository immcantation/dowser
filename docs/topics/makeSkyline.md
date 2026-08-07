**makeSkyline** - *get values for Bayesian Skyline plot*

Description
--------------------

`makeSkyline`


Usage
--------------------
```
makeSkyline(
object,
bins = 100,
youngest = 0,
clone_id = NULL,
max_height = c("min", "median", "mean", "max"),
exclude_germline = TRUE
)
```

Arguments
-------------------

object
:   treedata object with parameters_posterior and trees_posterior

bins
:   number of bins for plotting

youngest
:   timepoint of the most recently tip sampled (if 0, backward time used)

clone_id
:   name of the clone being analyzed (if desired)

max_height
:   max height to use (min, median, mean, max)

exclude_germline
:   exclude germline from skyline plot? (For TyCHE germline-rooted trees)




Value
-------------------

Bayesian Skyline values for given clone









