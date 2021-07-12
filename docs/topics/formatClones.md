**formatClones** - *Generate an ordered list of airrClone objects for lineage construction*

Description
--------------------

`formatClones` takes a `data.frame` or `tibble` with AIRR or 
Change-O style columns as input and masks gap positions, masks ragged ends, 
removes duplicates sequences, and merges annotations associated with duplicate
sequences. If specified, it will un-merge duplicate sequences with different 
values specified in the `trait` option. It returns a list of `airrClone`
objects ordered by number of sequences which serve as input for lineage reconstruction.


Usage
--------------------
```
formatClones(
data,
seq = "sequence_alignment",
clone = "clone_id",
subclone = "subclone_id",
nproc = 1,
chain = "H",
heavy = "IGH",
cell_id = "cell_id",
locus = "locus",
minseq = 2,
split_light = FALSE,
majoronly = FALSE,
columns = NULL,
...
)
```

Arguments
-------------------

data
:   data.frame containing the AIRR or Change-O data for a clone.
See [makeAirrClone](makeAirrClone.md) for required columns and their defaults

seq
:   sequence alignment column name.

clone
:   name of the column containing the identifier for the clone. All 
entries in this column should be identical.

subclone
:   name of the column containing the identifier for the subclone.

nproc
:   number of cores to parallelize formating over.

chain
:   if HL, include light chain information if available.

heavy
:   name of heavy chain locus (default = "IGH")

cell_id
:   name of the column containing cell assignment information

locus
:   name of the column containing locus information

minseq
:   minimum numbner of sequences per clone

split_light
:   split or lump subclones? See `getSubclones`.

majoronly
:   only return largest subclone and sequences without light chains

columns
:   additional data columns to include in output

...
:   additional arguments to pass to makeAirrClone




Value
-------------------

A tibble of [airrClone](airrClone-class.md) objects containing modified clones.


Details
-------------------

This function is a wrapper for [makeAirrClone](makeAirrClone.md). Also removes whitespace,
;, :, and = from ids



Examples
-------------------

```R
### Not run:
data(ExampleDb)
# clones <- formatClones(ExampleDb,trait="sample_id")
```



See also
-------------------

Executes in order [makeAirrClone](makeAirrClone.md). Returns a tibble of 
[airrClone](airrClone-class.md) objects 
which serve as input to [getTrees](getTrees.md) and [bootstrapTrees](bootstrapTrees.md).






