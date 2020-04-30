**ExampleDb** - *Example AIRR database*

Description
--------------------

A small example database subset from Laserson and Vigneault et al, 2014.


Usage
--------------------
```
ExampleDb
```




Format
-------------------

A data.frame with the following AIRR style columns:

+ `sequence_id`:           Sequence identifier
+ `sequence_alignment`:    IMGT-gapped observed sequence.
+ `germline_alignment_d_mask`:  IMGT-gapped germline sequence with N, P and 
D regions masked.
+ `v_call`:                V region allele assignments.
+ `v_call_genotyped`:      TIgGER corrected V region allele assignment.
+ `d_call`:                D region allele assignments.
+ `j_call`:                J region allele assignments.
+ `junction`:              Junction region sequence.
+ `junction_length`:       Length of the junction region in nucleotides.
+ `np1_length`:            Combined length of the N and P regions proximal
to the V region.
+ `np2_length`:            Combined length of the N and P regions proximal
to the J region.
+ `sample`:                Sample identifier. Time in relation to vaccination.
+ `isotype`:               Isotype assignment.
+ `duplicate_count`:       Copy count (number of duplicates) of the sequence.
+ `clone_id`:              Change-O assignment clonal group identifier.



References
-------------------


1. Laserson U and Vigneault F, et al. High-resolution antibody dynamics of 
vaccine-induced immune responses. 
Proc Natl Acad Sci USA. 2014 111:4928-33.





See also
-------------------

[ExampleDbChangeo](ExampleDbChangeo.md) [ExampleTrees](ExampleTrees.md)






