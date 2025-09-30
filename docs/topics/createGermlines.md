**createGermlines** - *[createGermlines](createGermlines.md) Determine consensus clone sequence and create germline for clone*

Description
--------------------

[createGermlines](createGermlines.md) Determine consensus clone sequence and create germline for clone


Usage
--------------------
```
createGermlines(
data,
references,
locus = "locus",
trim_lengths = FALSE,
force_trim = FALSE,
nproc = 1,
seq = "sequence_alignment",
v_call = "v_call",
d_call = "d_call",
j_call = "j_call",
amino_acid = FALSE,
id = "sequence_id",
clone = "clone_id",
v_germ_start = "v_germline_start",
v_germ_end = "v_germline_end",
v_germ_length = "v_germline_length",
d_germ_start = "d_germline_start",
d_germ_end = "d_germline_end",
d_germ_length = "d_germline_length",
j_germ_start = "j_germline_start",
j_germ_end = "j_germline_end",
j_germ_length = "j_germline_length",
np1_length = "np1_length",
np2_length = "np2_length",
na.rm = TRUE,
fields = NULL,
verbose = 0,
...
)
```

Arguments
-------------------

data
:   AIRR-table containing sequences from one clone

references
:   Full list of reference segments, see [readIMGT](readIMGT.md)

locus
:   Name of the locus column in the input data

trim_lengths
:   Remove trailing Ns from `seq` column if length different from germline?

force_trim
:   Remove all characters from sequence if different from germline? (not recommended)

nproc
:   Number of cores to use

seq
:   Column name for sequence alignment

v_call
:   Column name for V gene segment gene call

d_call
:   Column name for D gene segment gene call

j_call
:   Column name for J gene segment gene call

amino_acid
:   Perform reconstruction on amino acid sequence (experimental)

id
:   Column name for sequence ID

clone
:   Column name for clone ID

v_germ_start
:   Column name of index of V segment start within germline

v_germ_end
:   Column name of index of V segment end within germline

v_germ_length
:   Column name of index of V segment length within germline

d_germ_start
:   Column name of index of D segment start within germline

d_germ_end
:   Column name of index of D segment end within germline

d_germ_length
:   Column name of index of D segment length within germline

j_germ_start
:   Column name of index of J segment start within germline

j_germ_end
:   Column name of index of J segment end within germline

j_germ_length
:   Column name of index of J segment length within germline

np1_length
:   Column name in receptor specifying np1 segment length

np2_length
:   Column name in receptor specifying np2 segment length

na.rm
:   Remove clones with failed germline reconstruction?

fields
:   Character vector of additional columns to use for grouping. 
Sequences with disjoint values in the specified fields 
will be considered as separate clones.

verbose
:   amount of rubbish to print

...
:   Additional arguments passed to [buildGermline](buildGermline.md)




Value
-------------------

Tibble with reconstructed germlines


Details
-------------------

Return object adds/edits following columns:

+ `seq`:  Sequences potentially padded  same length as germline
+ `germline_alignment`: Full length germline
+ `germline_alignment_d_mask`: Full length, D region masked
+ `vonly`:   V gene segment of germline if vonly=TRUE
+ `regions`: String of VDJ segment in position if use_regions=TRUE




Examples
-------------------

```R
vdj_dir <- system.file("extdata", "germlines", "imgt", "human", "vdj", package="dowser")
imgt <- readIMGT(vdj_dir)

```


```
[1] "Read in 3 from 3 fasta files"

```


```R
db <- createGermlines(ExampleAirr[1,], imgt)
```

*Warning*:locus column not found, attempting to extract locus from V call*Warning*:Loci found: IGH*Warning*:Allele IGHV3-49*03 is not in the provided germline database.*Warning*:Allele IGHD6-19*01 is not in the provided germline database.*Warning*:Removing 1 failed clonal germlines. Clones: 3170

See also
-------------------

[createGermlines](createGermlines.md) [buildGermline](buildGermline.md), [stitchVDJ](stitchVDJ.md)






