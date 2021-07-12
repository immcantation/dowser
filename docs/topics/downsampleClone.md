**downsampleClone** - *`downsampleClone` Down-sample clone to maximum tip/switch ratio
TODO: Add support for weighting by collapseCount?*

Description
--------------------

`downsampleClone` Down-sample clone to maximum tip/switch ratio
TODO: Add support for weighting by collapseCount?


Usage
--------------------
```
downsampleClone(clone, trait, tip_switch = 20, tree = NULL)
```

Arguments
-------------------

clone
:   an [airrClone](airrClone-class.md) object

trait
:   trait considered for rarefaction
[getTrees](getTrees.md)

tip_switch
:   maximum tip/switch ratio

tree
:   a `phylo` tree object correspond to `clone`




Value
-------------------

A vector with sequence for each locus at a specified `node`
in `tree`.









