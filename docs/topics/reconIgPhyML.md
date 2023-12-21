**reconIgPhyML** - *Do IgPhyML maximum parsimony reconstruction*

Description
--------------------

`reconIgPhyML` IgPhyML parsimony reconstruction function


Usage
--------------------
```
reconIgPhyML(
file,
modelfile,
id,
igphyml = "igphyml",
mode = "switches",
type = "recon",
nproc = 1,
quiet = 0,
rm_files = FALSE,
rm_dir = NULL,
states = NULL,
palette = NULL,
resolve = 2,
rseed = NULL,
force_resolve = FALSE,
...
)
```

Arguments
-------------------

file
:   IgPhyML lineage file (see writeLineageFile)

modelfile
:   File specifying parsimony model

id
:   id for IgPhyML run

igphyml
:   location of igphyml executable

mode
:   return trees or count switches? (switches or trees)

type
:   get observed switches or permuted switches?

nproc
:   cores to use for parallelization

quiet
:   amount of rubbish to print

rm_files
:   remove temporary files?

rm_dir
:   remove temporary directory?

states
:   states in parsimony model

palette
:   deprecated

resolve
:   level of polytomy resolution. 0=none, 
1=maximum parsimony, 2=maximum ambiguity

rseed
:   random number seed if desired

force_resolve
:   continue even if polytomy resolution fails?

...
:   additional arguments




Value
-------------------

Either a tibble of switch counts or a list
of trees with internal nodes predicted by parsimony.









