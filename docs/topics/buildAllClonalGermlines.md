**buildAllClonalGermlines** - *[buildAllClonalGermlines](buildAllClonalGermlines.md) Determines and builds all possible germlines for a clone*

Description
--------------------

[buildAllClonalGermlines](buildAllClonalGermlines.md) Determines and builds all possible germlines for a clone


Usage
--------------------
```
buildAllClonalGermlines(
receptors,
references,
chain = "IGH",
use_regions = FALSE,
vonly = FALSE,
seq = "sequence_alignment",
id = "sequence_id",
clone = "clone_id",
v_call = "v_call",
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
j_germ_aa_length = "j_germline_aa_length",
amino_acid = FALSE,
threshold = 3,
...
)
```

Arguments
-------------------

receptors
:   AIRR-table containing sequences from one clone

references
:   Full list of reference segments, see [readIMGT](readIMGT.md)

chain
:   chain in `references` being analyzed

use_regions
:   Return string of VDJ regions? (optional)

vonly
:   Return germline of only v segment?

seq
:   Column name for sequence alignment

id
:   Column name for sequence ID

clone
:   Column name for clone ID

v_call
:   Column name for V gene segment gene call

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

j_germ_aa_length
:   Column name of J segment amino acid length (if amino_acid=TRUE)

amino_acid
:   Perform reconstruction on amino acid sequence (experimental)

threshold
:   The maximum germline length difference within a clone. Default is 3.

...
:   Additional arguments passed to [buildGermline](buildGermline.md)




Value
-------------------

A data frame with all possible reconstructed germlines




See also
-------------------

[createAllGermlines](createAllGermlines.md) [buildGermline](buildGermline.md), [stitchVDJ](stitchVDJ.md)






