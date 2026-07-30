**getTreesAndUCAs** - *[getTreesAndUCAs](getTreesAndUCAs.md) Construct trees and infer the UCA*

Description
--------------------

[getTreesAndUCAs](getTreesAndUCAs.md) Construct trees and infer the UCA


Usage
--------------------
```
getTreesAndUCAs(
clones,
data,
exec,
model_folder,
references,
dir = NULL,
model_folder_igk = NULL,
model_folder_igl = NULL,
partition = "single",
repertoire_wide = FALSE,
python = "python3",
id = "sample",
max_iters = 100,
nproc = 1,
rm_temp = TRUE,
quiet = 0,
chain = "H",
clone = "clone_id",
cell = "cell_id",
subsample_size = NA,
subsampling_method = c("random", "weighted", "least_mutated"),
search = c("codon", "nt"),
resolve_vj = FALSE,
fix_vj_in_cdr3 = TRUE,
fill_partials = TRUE,
split_light = FALSE,
...
)
```

Arguments
-------------------

clones
:   AIRR-table containing sequences [formatClones](formatClones.md)

data
:   The AIRR-table that was used to make the clones object.

exec
:   File path to the tree building executable

model_folder
:   The file path to the OLGA default model files for heavy chains

references
:   Reference genes. See [readIMGT](readIMGT.md)

dir
:   The file path of the directory of where data is saved. NULL is default.

model_folder_igk
:   The file path to the OLGA default model files for IGK

model_folder_igl
:   The file path to the OLGA default model files for IGL

partition
:   The partition model to use with IgPhyML. "single" is the default.

repertoire_wide
:   Build trees using parameters inferred from the entire dataset?

python
:   Specify the python call for your system. This is the call
on command line that issues the python you want to use. 
"python3" by default.

id
:   The run ID, sample by default

max_iters
:   The maximum number of iterations to run before ending. 
100 by default

nproc
:   The number of cores to use

rm_temp
:   Remove the generated files?

quiet
:   Amount of noise to print out

chain
:   Set to HL to use both heavy and light chain sequences

clone
:   The name of the clone id column used in [formatClones](formatClones.md).

cell
:   The name of the cell id in the AIRR table used to generate [formatClones](formatClones.md)

subsample_size
:   The amount that the clone should be sampled down to. 
By default this is NA to not induce subsampling.

subsampling_method
:   How to subsample. Methods include 'random', 'weighted', 
and 'least_mutated'. The later two methods
require 'mu_freq' to be passed as a trait when running

search
:   Search codon or nt space

resolve_vj
:   Resolve the V and J gene annotations within each clone?

fix_vj_in_cdr3
:   Check if the inferred V/J lengths go into the inferred cdr3 region and adjust accordingly.

fill_partials
:   A logical that will fill in the V and J UCAs of clones that have partial V/J sequence alignments

split_light
:   A logical that indicates if different light chain groups should be used to further split a clone (recommended for paired data)

...
:   Additional arguments passed to various other functions like [getTrees](getTrees.md) and [buildGermline](buildGermline.md)




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






