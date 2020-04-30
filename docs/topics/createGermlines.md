**createGermlines** - *Wrapper for CreateGermlines.py*

Description
--------------------

`createGermlines` reconstructs clonal and subclonal germlines


Usage
--------------------
```
createGermlines(
data,
exec,
refs,
file,
cf = "vj_clone",
format = "airr",
g = "dmask",
germ = "germline_alignment_d_mask",
rm_file = TRUE
)
```

Arguments
-------------------

data
:   tibble containing sequence information

exec
:   location of CreateGermline.py

refs
:   vector of reference allele locations

file
:   temporary file name to write

cf
:   column name for clone or subclone id

format
:   airr or changeo format

g
:   full or dmask germline option

germ
:   name of the column containg germline DNA sequences in
output file.

rm_file
:   remove temporary file?




Value
-------------------

a tibble containing reconstructed germlines









