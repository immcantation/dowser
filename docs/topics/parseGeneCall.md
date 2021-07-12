**parseGeneCall** - *Extract alleles from strings*

Description
--------------------

`parseGeneCall`


Usage
--------------------
```
parseGeneCall(gene, regex, action = "first")
```

Arguments
-------------------

gene
:   str containing V/D/J gene name

regex
:   regular expression string to extract allele

action
:   action to perform for multiple alleles;
one of ('first', 'set', 'list').




Value
-------------------

Allele name




See also
-------------------

[getVDJAllele](getVDJAllele.md)






