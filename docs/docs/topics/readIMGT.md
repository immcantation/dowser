**readIMGT** - *`readIMGT` read in IMGT database*

Description
--------------------

Loads all reference germlines from an Immcantation-formatted IMGT database.


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
Structure of object is list[[locus]][[segment]]
locus refers to locus (e.g. IGH, IGK, TRA)
segment refers to gene segment caegory (V, D, or J)


Details
-------------------

Input directory must be formatted to Immcantation standard.
See https://changeo.readthedocs.io/en/stable/examples/igblast.html for example
of how to download.



Examples
-------------------

```R
# vdj_dir contains a minimal example of reference germlines 
# (IGHV3-11*05, IGHD3-10*01 and IGHJ5*02)
# which are the gene assignments for ExamapleDb[1,]
vdj_dir <- system.file("extdata", "germlines", "imgt", "human", "vdj", package="dowser")
imgt <- readIMGT(vdj_dir)
```


```
[1] "Read in 3 from 3 fasta files"

```








