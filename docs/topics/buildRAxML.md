**buildRAxML** - *Wrapper to build RAxML-ng trees and infer intermediate nodes*

Description
--------------------

Wrapper to build RAxML-ng trees and infer intermediate nodes


Usage
--------------------
```
buildRAxML(
clone,
seq = "sequence",
exec,
model = "GTR",
partition = FALSE,
rseed = 28,
name = "run",
brln = "scaled",
starting_tree = NULL,
from_getTrees = FALSE,
rm_files = TRUE,
asr = TRUE,
rep = 1,
dir = NULL,
...
)
```

Arguments
-------------------

clone
:   list of `airrClone` objects

seq
:   the phylo_seq option does this clone uses. Possible options are "sequence", "hlsequence", or "lsequence"

exec
:   RAxML executable

model
:   The DNA model to be used. GTR is the default.

partition
:   A logical vector that indicates if you want to partition the data based on the heavy and light chains.

rseed
:   The random seed used for the parsimony inferences. This allows you to reproduce your results.

name
:   specifies the name of the output file

brln
:   A parameter that determines how branches are reported when partitioning. Options include scaled (default), 
unlinked, and linked

starting_tree
:   specifies a user starting tree file name and path in Newick format

from_getTrees
:   A logical that indicates if the desired starting tree is from getTrees and not a newick file

rm_files
:   remove temporary files?

asr
:   computes the marginal ancestral states of a tree

rep
:   Which repetition of the tree building is currently being run. Mainly for getBootstraps.

dir
:   Where the output files are to be made.

...
:   Additional arguments (not currently used)




Value
-------------------

`phylo` object created by RAxML-ng with nodes attribute
containing reconstructed sequences.









