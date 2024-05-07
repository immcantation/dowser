**maskSequences** - *`maskSequences` Mask codons split by insertions in V gene*

Description
--------------------

`maskSequences` Mask codons split by insertions in V gene


Usage
--------------------
```
maskSequences(
data,
sequence_id = "sequence_id",
sequence = "sequence",
sequence_alignment = "sequence_alignment",
v_sequence_start = "v_sequence_start",
v_sequence_end = "v_sequence_end",
v_germline_start = "v_germline_start",
v_germline_end = "v_germline_end",
junction_length = "junction_length",
keep_alignment = FALSE,
keep_insertions = FALSE,
mask_codons = TRUE,
mask_cdr3 = TRUE,
nproc = 1
)
```

Arguments
-------------------

data
:   BCR data table

sequence_id
:   sequence id column

sequence
:   input sequence column (query)

sequence_alignment
:   aligned (IMGT-gapped) sequence column (subject)

v_sequence_start
:   V gene start position in sequence

v_sequence_end
:   V gene end position in sequence

v_germline_start
:   V gene start position in sequence_alignment

v_germline_end
:   V gene end position in sequence_alignment

junction_length
:   name of junction_length column

keep_alignment
:   store alignment of query and subject sequences?

keep_insertions
:   return removed insertion sequences?

mask_codons
:   mask split codons?

mask_cdr3
:   mask CDR3 sequences?

nproc
:   number of cores to use




Value
-------------------

A tibble with masked sequence in sequence_masked column, 
 as well as other columns.


Details
-------------------

Performs global alignment of sequence and sequence_alignment, 
masking codons in sequence_alignment that are split by insertions (see examples)
masking_note notes codon positions in subject_alignment sequence that 
were masked, if found.
subject_alignment contains subject sequence aligned to query sequence (only 
if keep_alignment=TRUE)
query_alignment contains query sequence aligned to subject sequence (only if 
keep_alignment=TRUE)
sequence_masked will be NA if frameshift or alignment error detected. This 
will be noted
insertions column will be returned if keep_insertions=TRUE, contains a
comma-separated list of each <position in query alignment>-<sequence>. See example.
in masking_note.




See also
-------------------

[maskCodons](maskCodons.md), Biostrings::pairwiseAlignment.






