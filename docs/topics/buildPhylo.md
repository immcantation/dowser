**buildPhylo** - *Wrapper for alakazam::buildPhylipLineage*

Description
--------------------

Wrapper for alakazam::buildPhylipLineage


Usage
--------------------
```
buildPhylo(
clone,
exec,
temp_path = NULL,
verbose = 0,
rm_temp = TRUE,
seq = "sequence",
tree = NULL,
onetree = TRUE
)
```

Arguments
-------------------

clone
:   `airrClone` object

exec
:   dnapars or dnaml executable

temp_path
:   path to temporary directory

verbose
:   amount of rubbish to print

rm_temp
:   remove temporary files?

seq
:   sequence column in `airrClone` object

tree
:   fixed tree topology if desired (currently does nothing
if specified)

onetree
:   Only sample one tree if multiple found.




Value
-------------------

`phylo` object created by dnapars or dnaml with nodes attribute
containing reconstructed sequences.









