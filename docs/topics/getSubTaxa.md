**getSubTaxa** - *Get the tip labels as part of a clade defined by an internal node*

Description
--------------------

`getSubTaxa` Gets the tip labels from a clade


Usage
--------------------
```
getSubTaxa(node, tree)
```

Arguments
-------------------

node
:   node number that defines the target clade

tree
:   `phylo` object




Value
-------------------

A vector containing tip labels of the clade



Examples
-------------------

```R
# Get taxa from all subtrees
data(BiopsyTrees)
tree <- BiopsyTrees$trees[[8]]
all_subtrees <- lapply(1:length(tree$nodes), function(x)getSubTaxa(x, tree))

```








