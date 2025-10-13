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
iterations = 10,
ess_cutoff = 200,
ignore = c("traitfrequencies"),
quiet = 0,
...
)
```

Arguments
-------------------

clones
:   a tibble of `airrClone` objects, the output of
[formatClones](formatClones.md)

iterations
:   Maximum number of iterations

ess_cutoff
:   Minimum number of ESS for all parameters

ignore
:   Vector of parameters to ignore for ESS calculation

quiet
:   quiet notifications if > 0

...
:   Additional arguments for getTimeTrees




Value
-------------------

A tibble of `tidytree` and `airrClone` objects.


Details
-------------------

For examples and vignettes, see https://dowser.readthedocs.io









