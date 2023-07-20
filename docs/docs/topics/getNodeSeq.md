**getNodeSeq** - *Return IMGT gapped sequence of specified tree node*

Description
--------------------

`getNodeSeq` Sequence retrieval function.


Usage
--------------------
```
getNodeSeq(data, node, tree = NULL, clone = NULL, gaps = TRUE)
```

Arguments
-------------------

data
:   a tibble of `airrClone` objects, the output of 
[getTrees](getTrees.md)

node
:   numeric node in tree (see details)

tree
:   a `phylo` tree object containing `node`

clone
:   if `tree` not specified, supply clone ID in `data`

gaps
:   add IMGT gaps to output sequences?




Value
-------------------

A vector with sequence for each locus at a specified `node`
in `tree`.


Details
-------------------

Use plotTrees(trees)[[1]] + geom_label(aes(label=node))+geom_tippoint() to show
node labels, and getNodeSeq to return internal node sequences




See also
-------------------

[getTrees](getTrees.md)






