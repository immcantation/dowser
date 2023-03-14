**getSeq** - *Deprecated! Use getNodeSeq*

Description
--------------------

`getSeq` Sequence retrieval function.


Usage
--------------------
```
getSeq(data, node, tree = NULL, clone = NULL, gaps = TRUE)
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




See also
-------------------

[getTrees](getTrees.md)






