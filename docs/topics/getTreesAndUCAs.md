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
references = NULL,
clone = "clone_id",
cell = "cell_id",
heavy = "IGH",
sampling_method = "random",
subsample_size = NA,
search = "codon",
check_genes = TRUE,
igblast = NULL,
igblast_database = NULL,
ref_path = NULL,
organism = "human",
...
)
```

Arguments
-------------------

clones
:   AIRR-table containing sequences [formatClones](formatClones.md)

data
:   The AIRR-table that was used to make the clones object.

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

references
:   Reference genes. See [readIMGT](readIMGT.md)

clone
:   The name of the clone id column used in [formatClones](formatClones.md)

cell
:   The name of the cell id in the AIRR table used to generate [formatClones](formatClones.md)

heavy
:   The name of the heavy chain locus. Default is IGH.

sampling_method
:   How to subsample. Methods include 'random', 'lm' (least mutated),
and 'ratio' (a weighted sampling with heavier weights
towards the least mutated sequences). The later two methods
require 'mu_freq' to be passed as a trait when running
[formatClones](formatClones.md)

subsample_size
:   The amount that the clone should be sampled down to. Default is NA. Use NA if you do not wish to subsample -- in testing

search
:   Search codon or nt space

check_genes
:   Check if the inferred V/J lengths go into the inferred cdr3 region and adjust accordingly.

igblast
:   File path to igblast

igblast_database
:   The file path to the database setup for igblast

ref_path
:   The file path to your references parent folder.

organism
:   The type of organism to align to if using igblast.

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






