**resolvePolytomies** - *Resolve polytomies to have the minimum number of single timepoint clades*

Description
--------------------

Resolve polytomies to have the minimum number of single timepoint clades


Usage
--------------------
```
resolvePolytomies(
phy,
clone,
minlength = 0.001,
time = "time",
sequence = "sequence_id",
germline = "Germline",
verbose = FALSE
)
```

Arguments
-------------------

phy
:   Tree object

clone
:   airrClone data object corresponding to `phy`

minlength
:   Branch lengths to collapse in trees

time
:   Column name holding numeric time information

sequence
:   Column name holding sequence ID

germline
:   Germline sequence name

verbose
:   Print lots of rubbish while running?




Value
-------------------

A `phylo` tree object in which polytomies are resolved to 
have the minimum number of single timepoint clades.


Details
-------------------

Iteratively identifies polytomies (clusters of < minlength branches),
prunes each descendant branch, combines clades with the same timepoint
before grouping them back together. Checks to make sure that the divergence
of each tip is the same after resolution.




See also
-------------------

Uses output from [getTrees](getTrees.md) during [correlationTest](correlationTest.md).






