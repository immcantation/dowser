**calcRF** - *Finds the Robinson-Fould's cluster distance between phylogenies.*

Description
--------------------

`calcRF` Calculates the RF distance between two phylogenetic trees with 
the same tips and tip labels.


Usage
--------------------
```
calcRF(tree_1, tree_2, nproc = 1)
```

Arguments
-------------------

tree_1
:   A `phylo` object

tree_2
:   A `phylo` object

nproc
:   Number of cores to use for calculations.




Value
-------------------

The RF cluster value for the two input trees.









