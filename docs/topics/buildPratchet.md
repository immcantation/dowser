**buildPratchet** - *Wrapper for phangorn::pratchet*

Description
--------------------

Wrapper for phangorn::pratchet


Usage
--------------------
```
buildPratchet(
clone,
seq = "sequence",
asr = "seq",
asr_thresh = 0.05,
tree = NULL,
asr_type = "MPR",
verbose = 0,
resolve_random = TRUE,
data_type = "DNA"
)
```

Arguments
-------------------

clone
:   `airrClone` object

seq
:   sequence column in `airrClone` object

asr
:   return sequence or probability matrix?

asr_thresh
:   threshold for including a nucleotide as an alternative

tree
:   fixed tree topology if desired.

asr_type
:   MPR or ACCTRAN

verbose
:   amount of rubbish to print

resolve_random
:   randomly resolve polytomies?

data_type
:   Are sequences DNA or AA?




Value
-------------------

`phylo` object created by phangorn::pratchet with nodes
attribute containing reconstructed sequences.









