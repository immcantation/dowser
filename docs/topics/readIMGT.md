**readIMGT** - *[readIMGT](readIMGT.md) read in IMGT database
TODO: make auto-download or internal IMGT database*

Description
--------------------

[readIMGT](readIMGT.md) read in IMGT database
TODO: make auto-download or internal IMGT database


Usage
--------------------
```
readIMGT(dir, quiet = FALSE)
```

Arguments
-------------------

dir
:   directory containing Immcantation-formatted IMGT database

quiet
:   print warnings?




Value
-------------------

List of lists, leading to IMGT-gapped nucleotide sequences.
Structure of object is list[[organism]][[locus]][[segment]]
Organism refers to species (i.e. human, mouse)
locus refers to locus (e.g. IGH, IGK, TRA)
segment refers to gene segment caegory (V, D, or J)


Details
-------------------

Input directory must be formatted to Immcantation standard.
See https://changeo.readthedocs.io/en/stable/examples/igblast.html for example
of how to download.









