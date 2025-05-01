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
:   gap opening penalty (Biostrings::pairwiseAlignment)

gap_extension
:   gap extension penalty (Biostrings::pairwiseAlignment)

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

*Warning*:pairwiseAlignment() has moved to the pwalign package. Please call
  pwalign::pairwiseAlignment() to get rid of this warning.
```
$sequence_id
[1] 1

$sequence_masked
[1] ""

$masking_note
[1] ""

$insertions
[1] ""

$subject_alignment
[1] "ATC---ATCATC"

$query_alignment
[1] "ATCTTTATCATC"

$sequence_masked_v
[1] "ATCATCATC..."


```


```R

s <- "ATCATCATC..."
q <- "ATTTTCATCATC"
print(maskCodons("test",q,s,keep_alignment=TRUE,keep_insertions=TRUE))
```

*Warning*:pairwiseAlignment() has moved to the pwalign package. Please call
  pwalign::pairwiseAlignment() to get rid of this warning.
```
$sequence_id
[1] "test"

$sequence_masked
[1] ""

$masking_note
[1] "1,2"

$insertions
[1] "3-TTT"

$subject_alignment
[1] "AT---CATCATC"

$query_alignment
[1] "ATTTTCATCATC"

$sequence_masked_v
[1] "NNNATCATC..."


```



See also
-------------------

[maskSequences](maskSequences.md), Biostrings::pairwiseAlignment.






