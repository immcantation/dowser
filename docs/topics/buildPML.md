**buildPML** - *Wrapper for phangorn::optim.pml*

Description
--------------------

Wrapper for phangorn::optim.pml


Usage
--------------------
```
buildPML(
clone,
seq = "sequence",
sub_model = "GTR",
gamma = FALSE,
asr = "seq",
asr_thresh = 0.05,
tree = NULL,
data_type = "DNA",
optNni = TRUE,
optQ = TRUE,
optEdge = TRUE,
verbose = FALSE,
resolve_random = TRUE,
quiet = 0,
rep = NULL,
dir = NULL,
id = NULL,
asrp = FALSE
)
```

Arguments
-------------------

clone
:   `airrClone` object

seq
:   sequence column in `airrClone` object

sub_model
:   substitution model to use

gamma
:   gamma site rate variation?

asr
:   return sequence or probability matrix?

asr_thresh
:   threshold for including a nucleotide as an alternative

tree
:   fixed tree topology if desired.

data_type
:   Are sequences DNA or AA?

optNni
:   Optimize tree topology

optQ
:   Optimize Q matrix

optEdge
:   Optimize edge lengths

verbose
:   Print error messages as they happen?

resolve_random
:   randomly resolve polytomies?

quiet
:   amount of rubbish to print to console

rep
:   current bootstrap replicate (experimental)

dir
:   A directory to save the codon table

id
:   The identifier value

asrp
:   Get the codon table?




Value
-------------------

`phylo` object created by phangorn::optim.pml with nodes
attribute containing reconstructed sequences.









