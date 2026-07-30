**plotMSA** - *Simple function for plotting multiple sequence alignments
input can be either nucleotide or AA and contain ambiguous characters*

Description
--------------------

Simple function for plotting multiple sequence alignments
input can be either nucleotide or AA and contain ambiguous characters


Usage
--------------------
```
plotMSA(
data,
id,
seq,
seq_type = c("nt", "aa"),
text_size = 1.5,
palette = NULL
)
```

Arguments
-------------------

data
:   Data table with sequences

id
:   column name of sequence IDs (order will be preserved)

seq
:   column name of sequence data to be plotted

seq_type
:   sequence type in `seq` column, nucleotide (nt) or amino acid (aa)

text_size
:   size of sequence alignment text

palette
:   named vector of color for each character type in alignment, otherwise 
default nt or aa palette will be used




Value
-------------------

a ggplot2-formatted multiple sequence alignment plot




See also
-------------------

[getSeqPath](getSeqPath.md)






