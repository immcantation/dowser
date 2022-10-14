**buildIgphyml** - *Wrapper to build IgPhyML trees and infer intermediate nodes*

Description
--------------------

Wrapper to build IgPhyML trees and infer intermediate nodes


Usage
--------------------
```
buildIgphyml(
clone,
igphyml,
trees = NULL,
nproc = 1,
temp_path = NULL,
id = NULL,
rseed = NULL,
quiet = 0,
rm_files = TRUE,
rm_dir = NULL,
partition = c("single", "cf", "hl", "hlf", "hlc", "hlcf"),
omega = "e",
optimize = "lr",
motifs = "FCH",
hotness = "e,e,e,e,e,e",
asrc = 0.95,
splitfreqs = FALSE,
...
)
```

Arguments
-------------------

clone
:   list of `airrClone` objects

igphyml
:   igphyml executable

trees
:   list of tree topologies if desired

nproc
:   number of cores for parallelization

temp_path
:   path to temporary directory

id
:   IgPhyML run id

rseed
:   random number seed if desired

quiet
:   amount of rubbish to print

rm_files
:   remove temporary files?

rm_dir
:   remove temporary directory?

partition
:   How to partition omegas along sequences (see details)

omega
:   omega parameters to estimate (see IgPhyML docs)

optimize
:   optimize HLP rates (r), lengths (l), topology (t)

motifs
:   motifs to consider (see IgPhyML docs)

hotness
:   hotness parameters to estimate (see IgPhyML docs)

asrc
:   Intermediate sequence cutoff probability

splitfreqs
:   Calculate codon frequencies on each partition separately?

...
:   Additional arguments (not currently used)




Value
-------------------

`phylo` object created by igphyml with nodes attribute
containing reconstructed sequences.


Details
-------------------

Partition options:

+ `single`: 1 omega for whole sequence
+ `cf`: 2 omegas, 1 for all CDRs and 1 for all FWRs
+ `hl`: 2 omegas, 1 for heavy and 1 for light chain
+ `hlf`: 3 omegas, 1 for all CDRs, 2 for heavy/light FWRs
+ `hlc`: 3 omegas, 1 for all FWRs, 2 for heavy/light CDRs
+ `hlcf`: 4 omegas, 1 for each heavy/light CDR/FWR combination










