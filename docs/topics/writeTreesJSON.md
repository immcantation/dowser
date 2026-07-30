**writeTreesJSON** - *`writeTreesJSON`
Experimental. Write trees in AIRR format*

Description
--------------------

`writeTreesJSON`
Experimental. Write trees in AIRR format


Usage
--------------------
```
writeTreesJSON(
object,
file,
repertoire_id = "sample",
check = TRUE,
verbose = TRUE,
edge_tol = 1e-08,
cell = "cell_id",
heavy = "IGH",
dowser_fields = TRUE
)
```

Arguments
-------------------

object
:   Dowser object from getTrees

file
:   name of .json file

repertoire_id
:   repertoire_id to use

check
:   verify object is equivalent on reloading

verbose
:   print out more info

edge_tol
:   tolerance for branch length checks (if check=TRUE)

cell
:   cell id column name in Dowser object

heavy
:   name of heavy chain locus

dowser_fields
:   include dowser-specific information? (recommended)











