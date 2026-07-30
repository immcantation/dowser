**getDiffPoint** - *Recurse up to tree to find most recent node with different state, or the root*

Description
--------------------

`getDiffPoint`


Usage
--------------------
```
getDiffPoint(
tree,
targetnode,
trait,
height = "height",
verbose = FALSE,
eo_adjust = FALSE,
eo_type = NULL
)
```

Arguments
-------------------

tree
:   a treedata object, from getTimeTrees

targetnode
:   current node

trait
:   column name of the trait of interest

height
:   which height value to return

verbose
:   print out run info

eo_adjust
:   adjust heights using expectOccupancies (recommended, requires eo_type)

eo_type
:   if eo_adjust, trait value described by expectedOccupancies (typically state 1 of 2)











