**readTreesJSON** - *`readTreesJSON`
Experimental. Read trees from JSON/AIRR format from Dowser*

Description
--------------------

`readTreesJSON`
Experimental. Read trees from JSON/AIRR format from Dowser


Usage
--------------------
```
readTreesJSON(
file,
heavy = "IGH",
light = c("IGK", "IGL"),
verbose = TRUE,
edge_tol = 1e-08,
nproc = 1
)
```

Arguments
-------------------

file
:   .json file

heavy
:   name of heavy chain locus

light
:   names of light chain loci

verbose
:   how much info to print

edge_tol
:   tolerance for branch length checks (if check=TRUE)

nproc
:   number of cores to use (parallelizes by clone)




Details
-------------------

Reads files written by `[writeTreesJSON](writeTreesJSON.md)`, including trees built
with any `getTrees` `build` option.









