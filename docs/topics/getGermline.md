**getGermline** - *[getGermline](getGermline.md) get germline segment from specified receptor and segment*

Description
--------------------

[getGermline](getGermline.md) get germline segment from specified receptor and segment


Usage
--------------------
```
getGermline(
receptor,
references,
segment,
field,
germ_start,
germ_end,
germ_length,
germ_aa_start,
germ_aa_length,
amino_acid = FALSE
)
```

Arguments
-------------------

receptor
:   row from AIRR-table containing sequence of interest

references
:   list of reference segments. Must be specific to organism, 
locus, and segment

segment
:   Gene segment to search. Must be V, D, or J.

field
:   Column name for segment gene call (e.g. v_call)

germ_start
:   Column name of index of segment start within germline 
segment (e.g. v_germline_start)

germ_end
:   Similar to germ_start, but specifies end of segment 
(e.g. v_germline_end)

germ_length
:   Similar to germ_start, but specifies length of segment
(e.g. v_germline_end)

germ_aa_start
:   Column name of index of segment start within germline 
segment in AA (if amino_acid=TRUE, e.g. v_germline_start)

germ_aa_length
:   Similar to germ_start, but specifies length of segment
in AA (if amino_acid=TRUE, e.g. v_germline_end)

amino_acid
:   Perform reconstruction on amino acid sequence (experimental)




Value
-------------------

String of germline sequence from specified segment aligned with the 
sequence in the seq column of `receptor`.









