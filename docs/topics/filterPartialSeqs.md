**filterPartialSeqs** - *`filterPartialSeqs`*

Description
--------------------

`filterPartialSeqs`


Usage
--------------------
```
filterPartialSeqs(
data,
seq = "sequence_alignment",
cutoff = 250,
pattern = "[ATCG]",
verbose = TRUE
)
```

Arguments
-------------------

data
:   Data table of sequences (preferably in AIRR format)

seq
:   name of sequence column for filtering

cutoff
:   number of matches to the `pattern` option required

pattern
:   regex for matching (defaults to A, T, C, or G)

verbose
:   print out number of seqs removed
data(ExampleAirr)
filtered = filterPartialSeqs(ExampleAirr)











