**buildRAxML** - *Wrapper to build RAxML-ng trees and infer intermediate nodes*

Description
--------------------

Wrapper to build RAxML-ng trees and infer intermediate nodes


Usage
--------------------
```
buildRAxML(
clone,
exec,
seq = "sequence",
sub_model = "GTR",
partition = NULL,
rseed = 28,
name = "run",
starting_tree = NULL,
data_type = "DNA",
from_getTrees = FALSE,
rm_files = TRUE,
asr = TRUE,
rep = 1,
dir = NULL,
n_starts = NULL,
...
)
```

Arguments
-------------------

clone
:   list of `airrClone` objects

exec
:   RAxML-ng executable

seq
:   the phylo_seq option does this clone uses. Possible options are "sequence", "hlsequence", or "lsequence"

sub_model
:   The DNA model to be used. GTR is the default.

partition
:   A parameter that determines how branches are reported when partitioning. Options include NULL (default), 
scaled, unlinked, and linked

rseed
:   The random seed used for the parsimony inferences. This allows you to reproduce your results.

name
:   specifies the name of the output file

starting_tree
:   specifies a user starting tree file name and path in Newick format

data_type
:   Specifies what format your data is in, DNA or AA

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

n_starts
:   Number of max parsimony starting trees (default is 10 pars + 10 random)

...
:   Additional arguments (not currently used)




Value
-------------------

`phylo` object created by RAxML-ng with nodes attribute
containing reconstructed sequences.









