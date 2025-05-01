**getTreesAndUCAs** - *[getTreesAndUCAs](getTreesAndUCAs.md) Construct trees and infer the UCA*

Description
--------------------

[getTreesAndUCAs](getTreesAndUCAs.md) Construct trees and infer the UCA


Usage
--------------------
```
getTreesAndUCAs(
clones,
data = NULL,
dir = NULL,
build,
exec,
model_folder,
model_folder_igk = NULL,
model_folder_igl = NULL,
uca_script,
python = "python",
id = "sample",
max_iters = 100,
nproc = 1,
rm_temp = TRUE,
quiet = 0,
chain = "H",
omega = NULL,
optimize = "lr",
motifs = "FCH",
hotness = "e,e,e,e,e,e",
resolve_v = FALSE,
resolve_j = FALSE,
references = NULL,
clone = "clone_id",
...
)
```

Arguments
-------------------

clones
:   AIRR-table containing sequences [formatClones](formatClones.md)

data
:   The AIRR-table that was used to make the clones object. Required for resolve_v.

dir
:   The file path of the directory of where data is saved. NULL is default.

build
:   Name of the tree building method

exec
:   File path to the tree building executable

model_folder
:   The file path to the OLGA default model files for heavy chains

model_folder_igk
:   The file path to the OLGA default model files for IGK

model_folder_igl
:   The file path to the OLGA default model files for IGL

uca_script
:   The file path to the UCA python script

python
:   Specify the python call for your system. This is the call
on command line that issues the python you want to use. 
"python" by default.

id
:   The run ID

max_iters
:   The maximum number of iterations to run before ending

nproc
:   The number of cores to use

rm_temp
:   Remove the generated files?

quiet
:   Amount of noise to print out

chain
:   Set to HL to use both heavy and light chain sequences

omega
:   Omega parameters to estimate (see IgPhyMl docs)

optimize
:   Optimize HLP rates (r), lengths (l), and/or topology (r)

motifs
:   Motifs to consider (see IgPhyMl docs)

hotness
:   Hotness parameters to estimate (see IgPhyML docs)

resolve_v
:   Resolve the V gene as well

resolve_j
:   Resolve the J gene as well

references
:   Reference genes. See [readIMGT](readIMGT.md)

clone
:   The name of the clone id column used in [formatClones](formatClones.md)

...
:   Additional arguments passed to [buildGermline](buildGermline.md)




Value
-------------------

An `airrClone` object with trees and the inferred UCA


Details
-------------------

Return object adds/edits following columns:

+ `trees`:  The phylogenies associated with each clone
+ `UCA`:    The inferred UCA





See also
-------------------

[getTrees](getTrees.md)






