**stitchRegions** - *[stitchRegions](stitchRegions.md) Similar to [stitchVDJ](stitchVDJ.md) but with segment IDs 
instead of nucleotides*

Description
--------------------

[stitchRegions](stitchRegions.md) Similar to [stitchVDJ](stitchVDJ.md) but with segment IDs 
instead of nucleotides


Usage
--------------------
```
stitchRegions(
receptor,
v_seq,
d_seq,
j_seq,
np1_length = "np1_length",
np2_length = "np1_length",
n1_length = "n1_length",
p3v_length = "p3v_length",
p5d_length = "p5d_length",
p3d_length = "p3d_length",
n2_length = "n2_length",
p5j_length = "p5j_length",
np1_aa_length = "np1_aa_length",
np2_aa_length = "np2_aa_length",
amino_acid = FALSE
)
```

Arguments
-------------------

receptor
:   row from AIRR-table containing sequence of interest

v_seq
:   germline V segment sequence from [getGermline](getGermline.md)

d_seq
:   germline D segment sequence from [getGermline](getGermline.md)

j_seq
:   germline J segment sequence from [getGermline](getGermline.md)

np1_length
:   Column name in receptor specifying np1 segment length 
(e.g. np1_length)

np2_length
:   Column name in receptor specifying np2 segment length 
(e.g. np1_length)

n1_length
:   Column name in receptor specifying n1 segment length 
(experimental)

p3v_length
:   Column name in receptor specifying p3v segment length 
(experimental)

p5d_length
:   Column name in receptor specifying p5d segment length 
(experimental)

p3d_length
:   Column name in receptor specifying p3d segment length 
(experimental)

n2_length
:   Column name in receptor specifying n2 segment length 
(experimental)

p5j_length
:   Column name in receptor specifying p5j segment length 
(experimental)

np1_aa_length
:   Column name in receptor specifying np1 segment length 
in AA (if amino_acid=TRUE, e.g. np1_length)

np2_aa_length
:   Column name in receptor specifying np2 segment length 
in AA (if amino_acid=TRUE, e.g. np1_length)

amino_acid
:   Perform reconstruction on amino acid sequence (experimental)




Value
-------------------

Full length germline VDJ sequence with segment IDs instead of 
nucleotides.




See also
-------------------

[stitchVDJ](stitchVDJ.md)






