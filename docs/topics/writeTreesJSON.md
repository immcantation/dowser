**writeTreesJSON** - *`writeTreesJSON`
Experimental. Write trees in AIRR format*

Description
--------------------

`writeTreesJSON`
Experimental. Write trees in AIRR format


Usage
--------------------
```
writeTreesJSON(
object,
file,
repertoire_id = "sample",
check = TRUE,
verbose = TRUE,
edge_tol = 1e-08,
cell = "cell_id",
heavy = "IGH",
light = c("IGK", "IGL"),
dowser_fields = TRUE,
nproc = 1
)
```

Arguments
-------------------

object
:   Dowser object from getTrees

file
:   name of .json file

repertoire_id
:   repertoire_id to use

check
:   verify object is equivalent on reloading

verbose
:   print out more info

edge_tol
:   tolerance for branch length checks (if check=TRUE)

cell
:   cell id column name in Dowser object

heavy
:   name of heavy chain locus

light
:   names of light chain loci

dowser_fields
:   include dowser-specific information? (recommended)

nproc
:   number of cores to use (parallelizes by clone)




Details
-------------------

Works with trees built by any of `getTrees`'s `build` options
(`"pratchet"`, `"pml"`, `"igphyml"`, `"raxml"`).
`getTrees(..., build="pml")` trees store the full
`phangorn::optim.pml` fit in `tree$parameters`; since that object
isn't JSON-serializable (and isn't meaningful to reconstruct from a file),
it's reduced to a flat list of fitted model parameters via
`pmlParamsToList` before being written out. See that function's
documentation for what is kept and why.









