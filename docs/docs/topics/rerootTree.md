**rerootTree** - *Reroot phylogenetic tree to have its germline sequence at a zero-length branch
to a node which is the direct ancestor of the tree's UCA. Assigns `uca`
to be the ancestral node to the tree's germline sequence, as `germid` as
the tree's germline sequence ID.*

Description
--------------------

Reroot phylogenetic tree to have its germline sequence at a zero-length branch
to a node which is the direct ancestor of the tree's UCA. Assigns `uca`
to be the ancestral node to the tree's germline sequence, as `germid` as
the tree's germline sequence ID.


Usage
--------------------
```
rerootTree(tree, germline, min = 0.001, verbose = 1)
```

Arguments
-------------------

tree
:   An ape `phylo` object

germline
:   ID of the tree's predicted germline sequence

min
:   Maximum allowed branch length from germline to root

verbose
:   amount of rubbish to print




Value
-------------------

`phylo` object rooted at the specified germline









