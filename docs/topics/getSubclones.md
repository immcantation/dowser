**getSubclones** - *#' Deprecated! Use resolveLightChains*

Description
--------------------

`getSubClones` plots a tree or group of trees


Usage
--------------------
```
getSubclones(
heavy,
light,
nproc = 1,
minseq = 1,
id = "sequence_id",
seq = "sequence_alignment",
clone = "clone_id",
cell = "cell_id",
v_call = "v_call",
j_call = "j_call",
junc_len = "junction_length",
nolight = "missing"
)
```

Arguments
-------------------

heavy
:   a tibble containing heavy chain sequences with clone_id

light
:   a tibble containing light chain sequences

nproc
:   number of cores for parallelization

minseq
:   minimum number of sequences per clone

id
:   name of the column containing sequence identifiers.

seq
:   name of the column containing observed DNA sequences. All 
sequences in this column must be multiple aligned.

clone
:   name of the column containing the identifier for the clone. All 
entries in this column should be identical.

cell
:   name of the column containing identifier for cells.

v_call
:   name of the column containing V-segment allele assignments. All 
entries in this column should be identical to the gene level.

j_call
:   name of the column containing J-segment allele assignments. All 
entries in this column should be identical to the gene level.

junc_len
:   name of the column containing the length of the junction as a 
numeric value. All entries in this column should be identical 
for any given clone.

nolight
:   string to use to indicate a missing light chain




Value
-------------------

a tibble containing









