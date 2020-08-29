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
...
)
```

Arguments
-------------------

clone
:   `airrClone` object

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

...
:   Additional arguments (not currently used)




Value
-------------------

`phylo` object created by igphyml with nodes attribute
containing reconstructed sequences.









