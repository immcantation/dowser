**checkNodeDivergences** - *For each node in nodes list, check that the current divergence
is the same as previously reported*

Description
--------------------

For each node in nodes list, check that the current divergence
is the same as previously reported


Usage
--------------------
```
checkNodeDivergences(tree, stop = TRUE, catch_null = TRUE)
```

Arguments
-------------------

tree
:   A `phylo` object with nodes list

stop
:   Throw an error if divergences not equal

catch_null
:   Should null divergence values count as mismatches?




Value
-------------------

TRUE or FALSE depening on whether node divergences are consistent
if stop=TRUE, a FALSE value will cause an error









