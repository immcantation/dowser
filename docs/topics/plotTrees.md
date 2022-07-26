**plotTrees** - *Plot a tree with colored internal node labels using ggtree*

Description
--------------------

`plotTrees` plots a tree or group of trees


Usage
--------------------
```
plotTrees(
trees,
nodes = FALSE,
tips = NULL,
tipsize = NULL,
scale = 0.01,
node_palette = "Dark2",
tip_palette = node_palette,
base = FALSE,
layout = "rectangular",
node_nums = FALSE,
tip_nums = FALSE,
title = TRUE,
labelsize = NULL,
common_scale = FALSE,
ambig = "blend"
)
```

Arguments
-------------------

trees
:   A tibble containing `phylo` and `airrClone`
objects

nodes
:   color internal nodes if possible?

tips
:   color tips if possible?

tipsize
:   size of tip shape objects

scale
:   width of branch length scale bar

node_palette
:   color palette for nodes

tip_palette
:   color palette for tips

base
:   recursion base case (don't edit)

layout
:   rectangular or circular tree layout?

node_nums
:   plot internal node numbers?

tip_nums
:   plot tip numbers?

title
:   use clone id as title?

labelsize
:   text size

common_scale
:   strecth plots so branches are on same scale?
determined by sequence with highest divergence

ambig
:   How to color ambiguous node reconstructions? (blend or grey)




Value
-------------------

a grob containing a tree plotted by `ggtree`.


Details
-------------------

Function uses `ggtree` functions to plot tree topologlies estimated by 
[getTrees](getTrees.md), and [findSwitches](findSwitches.md). Object can be further modified with 
`ggtree` functions. Please check out 
https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/ggtree.html and
cite `ggtree` in addition to `dowser` if you use this function.



Examples
-------------------

```R
data(ExampleClones)
trees <- getTrees(ExampleClones[10,])
plotTrees(trees)[[1]]
```

![2](plotTrees-2.png)


See also
-------------------

[getTrees](getTrees.md), [findSwitches](findSwitches.md)






