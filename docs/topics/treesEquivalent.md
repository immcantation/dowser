**treesEquivalent** - *Check whether two tree objects are equivalent*

Description
--------------------

Check whether two tree objects are equivalent


Usage
--------------------
```
treesEquivalent(
obja,
objb,
edge_tol = 1e-08,
numbering_match = FALSE,
clonesa = NULL,
clonesb = NULL,
check_extended = FALSE,
gaps = TRUE
)
```

Arguments
-------------------

obja
:   First phylo or treedata object

objb
:   Second phylo or treedata object

edge_tol
:   tolerance for branch length checks (if check=TRUE)

numbering_match
:   require internal node numbers to match?

clonesa
:   Dowser clones object associated with obja (if check_extended=TRUE)

clonesb
:   Dowser clones object associated with objb (if check_extended=TRUE)

check_extended
:   Also check node sequences and state vector? Requires clonesa/b to be supplied

gaps
:   check sequences wth IMGT gaps with check_extended?




Details
-------------------

For treedata objects, check both @phylo and @data









