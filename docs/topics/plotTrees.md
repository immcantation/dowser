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
labelsize = NULL
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




Value
-------------------

a grob containing a tree plotted by `ggtree`.


Details
-------------------

Function uses `ggtree` functions to plot tree topologlies estimated by 
[getTrees](getTrees.md), and [bootstrapTrees](bootstrapTrees.md). Object can be further modified with 
`ggtree` functions. Please check out 
https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/ggtree.html and
cite `ggtree` in addition to `dowser` if you use this function.



Examples
-------------------

```R
### Not run:
data(ExampleAirr)
# ExampleAirr$sample_id <- sample(ExampleAirr$sample_id)
# clones <- formatClones(ExampleAirr, trait="sample_id")
# 
# trees <- getTrees(clones[1:2,])
# plots <- plotTrees(trees)
# plots[[1]]
```



See also
-------------------

[getTrees](getTrees.md), [bootstrapTrees](bootstrapTrees.md)






