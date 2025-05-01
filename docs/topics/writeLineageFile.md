**writeLineageFile** - *Write lineage file for IgPhyML use*

Description
--------------------

Write lineage file for IgPhyML use


Usage
--------------------
```
writeLineageFile(
data,
trees = NULL,
dir = ".",
id = "N",
rep = NULL,
trait = NULL,
empty = TRUE,
partition = "single",
heavy = "IGH",
...
)
```

Arguments
-------------------

data
:   list of `airrClone` objects

trees
:   list of `phylo` objects corresponding to `data`

dir
:   directory to write file

id
:   id used for IgPhyML run

rep
:   bootstrap replicate

trait
:   string appended to sequence id in fasta files

empty
:   output uninformative sequences?

partition
:   how to partition omegas

heavy
:   name of heavy chain locus

...
:   additional arguments to be passed




Value
-------------------

Name of created lineage file.









