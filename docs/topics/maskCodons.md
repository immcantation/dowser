**maskCodons** - *`maskCodons` Masks codons split by insertions*

Description
--------------------

`maskCodons` Masks codons split by insertions


Usage
--------------------
```
maskCodons(
id,
q,
s,
keep_alignment = FALSE,
gap_opening = 5,
gap_extension = 1,
keep_insertions = FALSE,
mask = TRUE
)
```

Arguments
-------------------

id
:   sequence id

q
:   (query) un-aligned input sequence (sequence)

s
:   (subject) aligned input sequence (sequence_alignment)

keep_alignment
:   store q and s alignments

gap_opening
:   gap opening penalty (pwalign::pairwiseAlignment)

gap_extension
:   gap extension penalty (pwalign::pairwiseAlignment)

keep_insertions
:   return removed insertion sequences?

mask
:   if FALSE, don't mask codons




Value
-------------------

A list with split codons masked, if found (sequence_masked).


Details
-------------------

Performs global alignment of q and s, masks codons in s that are split by 
insertions (see example)
masking_note notes codon positions in subject_alignment sequence that were 
masked, if found.
subject_alignment contains subject sequence aligned to query (q) sequence
query_alignment contains query sequence aligned to subject (q) sequence
sequence_masked will be NA if frameshift or alignment error detected/



Examples
-------------------

```R
s = "ATCATCATC..."
q = "ATCTTTATCATC"
print(maskCodons(1,q,s,TRUE))

```

**Error in loadNamespace(x)**: there is no package called ‘pwalign’
```R

s <- "ATCATCATC..."
q <- "ATTTTCATCATC"
print(maskCodons("test",q,s,keep_alignment=TRUE,keep_insertions=TRUE))

```

**Error in loadNamespace(x)**: there is no package called ‘pwalign’

See also
-------------------

[maskSequences](maskSequences.md), pwalign::pairwiseAlignment.






