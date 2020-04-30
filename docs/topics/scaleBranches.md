**scaleBranches** - *Scale branch lengths to represent either mutations or mutations per site.*

Description
--------------------

`scaleBranches` Branch length scaling function.


Usage
--------------------
```
scaleBranches(clones, edge_type = "mutations")
```

Arguments
-------------------

clones
:   a tibble of `airrClone` and `phylo` objects,
the output of [getTrees](getTrees.md).

edge_type
:   Either `genetic_distance` (mutations per site) or 
`mutations`




Value
-------------------

A tibble with `phylo` objects that have had branch lengths 
rescaled as specified.


Details
-------------------

Uses clones$trees[[1]]$edge_type to determine how branches are currently scaled.




See also
-------------------

[getTrees](getTrees.md)






