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
model = "GTR",
gamma = FALSE,
asr = "seq",
asr_thresh = 0.05,
tree = NULL,
...
)
```

Arguments
-------------------

clone
:   `airrClone` object

seq
:   sequece column in `airrClone` object

model
:   substitution model to use

gamma
:   gamma site rate variation?

asr
:   return sequence or probability matrix?

asr_thresh
:   threshold for including a nucleotide as an alternative

tree
:   fixed tree topology if desired.

...
:   Additional arguments (not currently used)




Value
-------------------

`phylo` object created by phangorn::optim.pml with nodes
attribute containing reconstructed sequences.









