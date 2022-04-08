**getDivergence** - *Get divergence from root of tree for each tip*

Description
--------------------

`getDivergence` get sum of branch lengths leading from the 
root of the tree. If the germline sequence is included in the tree,
this will equal the germline divergence. If germline removed,
this will equal the MRCA divergence


Usage
--------------------
```
getDivergence(phy, minlength = 0.001)
```

Arguments
-------------------

phy
:   Tree object

minlength
:   Branch lengths to collapse in trees




Value
-------------------

A named vector of each tip's divergence from the tree's root.









