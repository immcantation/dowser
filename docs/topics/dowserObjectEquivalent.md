**dowserObjectEquivalent** - *`dowserObjectEquivalent`
Experimental. Check if two Dowser objects are equivalent*

Description
--------------------

`dowserObjectEquivalent`
Experimental. Check if two Dowser objects are equivalent


Usage
--------------------
```
dowserObjectEquivalent(
obj1,
obj2,
verbose = TRUE,
edge_tol = 1e-08,
dowser_fields = TRUE,
nproc = 1
)
```

Arguments
-------------------

obj1
:   First Dowser object

obj2
:   Second Dowser object

verbose
:   print out more info

edge_tol
:   tolerance for branch length checks (if check=TRUE)

dowser_fields
:   check dowser-specific fields and gapped sequences?

nproc
:   number of cores to use




Details
-------------------

In addition to the existing tree topology, edge length, sequence, and
data slot checks, this also verifies `tree$parameters` when present
-- including `build="pml"` trees.









