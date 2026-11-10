**refitUCA** - *[refitUCA](refitUCA.md) Refit the tree and ancestral sequence reconstruction using the inferred UCA*

Description
--------------------

[refitUCA](refitUCA.md) Refit the tree and ancestral sequence reconstruction using the inferred UCA


Usage
--------------------
```
refitUCA(
clones,
exec,
dir = NULL,
partition = "single",
id = "sample",
nproc = 1,
rm_temp = TRUE,
quiet = 0,
...
)
```

Arguments
-------------------

clones
:   The output of [getTreesAndUCAs](getTreesAndUCAs.md)

exec
:   File path to the IgPhyML executable

dir
:   The file path of the directory of where data is saved. NULL is default.

partition
:   The partition model to use with IgPhyML. "single" is the default.

id
:   The run ID, sample by default

nproc
:   The number of cores to use

rm_temp
:   Remove the generated files?

quiet
:   Amount of noise to print out

...
:   Additional arguments passed to [getTrees](getTrees.md)




Value
-------------------

An `airrClone` object with trees refit using the UCA as the germline


Details
-------------------

The UCA stored in the germline slots by [getTreesAndUCAs](getTreesAndUCAs.md) is used
as the germline for each clone and the trees are rebuilt with IgPhyML. The
germline slots hold the UCA trimmed to the sites in the tip sequences, rather
than the full length UCA in `UCA$ungapped`. The germline node of each refit
tree is set to this UCA and the UCA log-likelihoods are carried over to the new
tree parameters. The `UCA` column is not changed.




See also
-------------------

[getTreesAndUCAs](getTreesAndUCAs.md) [getTrees](getTrees.md)






