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
omega = NULL,
optimize = "lr",
motifs = "FCH",
hotness = "e,e,e,e,e,e",
rates = NULL,
asrc = 0.95,
splitfreqs = FALSE,
asrp = FALSE,
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

rates
:   comma delimited list showing which omega-defined partitions
get a separate rate (e.g. omega=e,e rates=0,1).

asrc
:   Intermediate sequence cutoff probability

splitfreqs
:   Calculate codon frequencies on each partition separately?

asrp
:   Run ASRp?

...
:   Additional arguments (not currently used)




Value
-------------------

`phylo` object created by igphyml with nodes attribute
containing reconstructed sequences.


Details
-------------------

Partition options in rate order:

+ `single`: 1 omega for whole sequence
+ `cf`: 2 omegas, 1 for all CDRs and 1 for all FWRs
+ `hl`: 2 omegas, 1 for heavy and 1 for light chain
+ `hlf`: 3 omegas, 1 for heavy FWR, 1 for all CDRs, and 1 for light FWRs
+ `hlc`: 3 omegas, 1 for all FWRs, 1 for heavy CDRs, and 1 for light CDRs
+ `hlcf`: 4 omegas, 1 for each heavy FWR, 1 for heavy CDR, 1 for light FWR, and 1 for light CDR










