**colorTrees** - *Get a color palette for a predefined set of trait values*

Description
--------------------

`colorTree` Gets a color palette for a predefined set of trait values


Usage
--------------------
```
colorTrees(trees, palette, ambig = "blend")
```

Arguments
-------------------

trees
:   list of phylo objects with assigned internal node states

palette
:   named vector of colors (see [getPalette](getPalette.md))

ambig
:   how should ambiguous states be colored (blend or grey)




Value
-------------------

A list of colored trees


Details
-------------------

Trees must have node states represented in a "states" vector. By default,
ambiguous states (separated by ",") have their colors blended. If




See also
-------------------

[getPalette](getPalette.md), [getTrees](getTrees.md), [plotTrees](plotTrees.md)






