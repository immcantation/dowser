**buildGermline** - *`buildGermline` reconstruct germline segments from alignment data*

Description
--------------------

Reconstruct germlines from alignment data.


Usage
--------------------
```
buildGermline(
receptor,
references,
seq = "sequence_alignment",
id = "sequence_id",
clone = "clone_id",
v_call = "v_call",
d_call = "d_call",
j_call = "j_call",
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
amino_acid = FALSE
)
```

Arguments
-------------------

receptor
:   row from AIRR-table containing sequence of interest

references
:   list of reference segments. Must be specific to organism 
and locus

seq
:   Column name for sequence alignment

id
:   Column name for sequence ID

clone
:   Column name for clone ID

v_call
:   Column name for V gene segment gene call

d_call
:   Column name for D gene segment gene call

j_call
:   Column name for J gene segment gene call

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

amino_acid
:   Perform reconstruction on amino acid sequence (experimental)




Value
-------------------

List of reconstructed germlines


Details
-------------------

Return object contains multiple IMGT-gapped germlines:

+ `full`:    Full length germline
+ `dmask`:   Full length germline with D region masked
+ `vonly`:   V gene segment of germline
+ `regions`: String showing VDJ segment of each position





See also
-------------------

[buildClonalGermline](buildClonalGermline.md), [stitchVDJ](stitchVDJ.md)






