**buildClonalGermline** - *`buildClonalGermline` Determine consensus clone sequence and create germline for clone*

Description
--------------------

Determine consensus clone sequence and create germline for clone


Usage
--------------------
```
buildClonalGermline(
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
j_germ_length = "j_germline_length",
j_germ_aa_length = "j_germline_aa_length",
amino_acid = FALSE,
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

j_germ_length
:   Column name of J segment length within germline

j_germ_aa_length
:   Column name of J segment amino acid length (if amino_acid=TRUE)

amino_acid
:   Perform reconstruction on amino acid sequence (experimental)

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





See also
-------------------

[createGermlines](createGermlines.md) [buildGermline](buildGermline.md), [stitchVDJ](stitchVDJ.md)






