**dfToFasta** - *Write a fasta file of sequences
`readFasta` reads a fasta file*

Description
--------------------

Write a fasta file of sequences
`readFasta` reads a fasta file


Usage
--------------------
```
dfToFasta(
df,
file,
id = "sequence_id",
seq = "sequence",
imgt_gaps = FALSE,
columns = NULL
)
```

Arguments
-------------------

df
:   dataframe of sequences

file
:   FASTA file for output

id
:   Column name of sequence ids

seq
:   Column name of sequences

imgt_gaps
:   Keep IMGT gaps if present?

columns
:   vector of column names to append to sequence id




Value
-------------------

File of FASTA formatted sequences









