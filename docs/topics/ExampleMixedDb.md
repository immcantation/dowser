**ExampleMixedDb** - *Example Change-O database*

Description
--------------------

A small example database subset from Turner, J. S. et al. Human germinal centres 
engage memory and naive B cells after influenza vaccination. Nature 586, 127–132 (2020).


Usage
--------------------
```
ExampleMixedDb
```




Format
-------------------

A data.frame with the following Change-O style columns:

+ `sequence_id`:           Sequence identifier
+  `sequence`:               B cell sequence
+  `productive`:             A logical indicating if the sequence is productive.
+ `v_call`:                V region allele assignments.
+  `d_call`:                 D region allele assignments. 
+ `j_call`:                J region allele assignments.
+ `sequence_alignment`:    Sequence alignment.
+  `germline_alignment`:     Germline alignment without gaps. 
+  `junction`:               Junction
+  `juncation_aa`:           Junction aa
+  `vj_inframe`:             A logical to see if the vj genes are in frame
+  `stop_codon`:             A indicator if there is a stop codon within the alignment
+  `locus`:                  Locus identifier. 
+  `v_sequence_start`:       Where the V gene starts
+  `v_sequence_end`:         Where the V gene ends
+  `v_germline_start`:       Where the V germline starts
+  `v_germline_end`:         Where the V germline ends
+  `np1_length`:             Length of np1
+  `d_sequence_start`:       Where the D gene starts
+  `d_sequence_end`:         Where the D gene ends
+  `d_germline_start`:       Where the D germline starts
+  `d_germline_end`:         Where the D germline ends
+  `np2_length`:             Length of np2
+  `j_sequence_start`:       Where the J gene starts
+  `j_sequence_end`:         Where the J gene ends
+  `j_germline_start`:       Where the J germline starts
+  `j_germline_end`:         Where the J germline ends
+ `junction_length`:       Length of the junction region in nucleotides.
+  `v_score`:                V score
+  `v_identity`:             Identity score of V
+  `v_support`:              V support
+  `d_score`:                D score
+  `d_identity`:             D identity 
+  `d_support`:              D support
+  `j_score`:                J score
+  `j_support`:              J support
+  `j_identity`:             J identity 
+  `cell_id`:                Cell identifier 
+  `consensus_count`:        Consensus count 
+  `indels`:                 Logical if indels are present 
+  `sequence_vdj`:           VDJ sequence
+  `v_germ_start_vdj`:       Where the V germline starts on the VDJ
+  `v_germ_end_vdj`:         Where the V germline ends on the VDJ
+  `subject`:                Subject identifier 
+  `timepoint`:              Day the sample was taken 
+  `cell_type`:              Type of cell 
+  `replicate`:              Replicate number 
+ `clone_id`:              Change-O assignment clonal group identifier.
+  `seq_type`:               Identifier of data type (10x)
+  `vj_gene`:                VJ gene
+  `vj_alt_gene`:            Alternative VJ gene
+  `v_germline_length`:      Length of the V germline segment
+  `d_germline_length`:      Length of the D germline segment 
+  `j_germline_lenght`:      Length of the J germline segment 
+  `germline_alignment_d_mask`:  Germline alignment with gaps










