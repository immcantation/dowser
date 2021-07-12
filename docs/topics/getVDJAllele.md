**getVDJAllele** - *`getVDJAllele` Extract allele from gene call string*

Description
--------------------

`getVDJAllele` Extract allele from gene call string


Usage
--------------------
```
getVDJAllele(gene, segment, action = "first")
```

Arguments
-------------------

gene
:   string with V gene calls

segment
:   Gene segment to search. Must be V, D, or J.

action
:   action to perform for multiple alleles;
one of ('first', 'set', 'list').




Value
-------------------

Allele call for given gene segment




See also
-------------------

[parseGeneCall](parseGeneCall.md), [buildClonalGermline](buildClonalGermline.md)






