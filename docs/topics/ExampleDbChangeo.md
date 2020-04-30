**ExampleDbChangeo** - *Example Change-O database*

Description
--------------------

A small example database subset from Laserson and Vigneault et al, 2014.


Usage
--------------------
```
ExampleDbChangeo
```




Format
-------------------

A data.frame with the following Change-O style columns:

+ `SEQUENCE_ID`:           Sequence identifier
+ `SEQUENCE_IMGT`:         IMGT-gapped observed sequence.
+ `GERMLINE_IMGT_D_MASK`:  IMGT-gapped germline sequence with N, P and 
D regions masked.
+ `V_CALL`:                V region allele assignments.
+ `V_CALL_GENOTYPED`:      TIgGER corrected V region allele assignment.
+ `D_CALL`:                D region allele assignments.
+ `J_CALL`:                J region allele assignments.
+ `JUNCTION`:              Junction region sequence.
+ `JUNCTION_LENGTH`:       Length of the junction region in nucleotides.
+ `NP1_LENGTH`:            Combined length of the N and P regions proximal
to the V region.
+ `NP2_LENGTH`:            Combined length of the N and P regions proximal
to the J region.
+ `SAMPLE`:                Sample identifier. Time in relation to vaccination.
+ `ISOTYPE`:               Isotype assignment.
+ `DUPCOUNT`:              Copy count (number of duplicates) of the sequence.
+ `CLONE`:                 Change-O assignment clonal group identifier.



References
-------------------


1. Laserson U and Vigneault F, et al. High-resolution antibody dynamics of 
vaccine-induced immune responses. 
Proc Natl Acad Sci USA. 2014 111:4928-33.





See also
-------------------

[ExampleDb](ExampleDb.md) [ExampleTrees](ExampleTrees.md)






