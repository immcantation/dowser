**collapseNodes** - *Collapse internal nodes with the same predicted sequence*

Description
--------------------

`collapseNodes` Node collapsing function.


Usage
--------------------
```
collapseNodes(trees, tips = FALSE, check = TRUE)
```

Arguments
-------------------

trees
:   a tibble of `airrClone` objects, the output of [getTrees](getTrees.md)

tips
:   collapse tips to internal nodes? (experimental)

check
:   check that collapsed nodes are consistent with original tree




Value
-------------------

A tibble with `phylo` objects that have had internal nodes collapsed.


Details
-------------------

Use plotTrees(trees)[[1]] + geom_label(aes(label=node)) + geom_tippoint() to show
node labels, and getSeq to return internal node sequences




See also
-------------------

[getTrees](getTrees.md)






