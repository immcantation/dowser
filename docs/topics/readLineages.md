**readLineages** - *Read in all trees from a lineages file*

Description
--------------------

Read in all trees from a lineages file


Usage
--------------------
```
readLineages(
file,
states = NULL,
palette = "Dark2",
run_id = "",
quiet = TRUE,
append = NULL,
format = "nexus",
type = "jointpars"
)
```

Arguments
-------------------

file
:   IgPhyML lineage file

states
:   states in parsimony model

palette
:   palette for coloring internal nodes

run_id
:   id used for IgPhyML run

quiet
:   avoid printing rubbish on screen?

append
:   string appended to fasta files

format
:   format of input file with trees

type
:   Read in parsimony reconstructions or ancestral sequence
reconstructions? "jointpars" reads in parsimony states, 
others read in sequences in internal nodes




Value
-------------------

A list of phylo objects from `file`.









