**stitchVDJ** - *[stitchVDJ](stitchVDJ.md) combines germline gene segments to a single string*

Description
--------------------

[stitchVDJ](stitchVDJ.md) combines germline gene segments to a single string


Usage
--------------------
```
stitchVDJ(
receptor,
v_seq,
d_seq,
j_seq,
np1_length = "np1_length",
np2_length = "np2_length",
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

Full length germline VDJ sequence aligned with aligned with the 
sequence in the `seq` column of `receptor`.









