**getSubclones** - *Define subclones based on light chain rearrangements*

Description
--------------------

`getSubclones` plots a tree or group of trees


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
cell_id = "cell_id",
v_call = "v_call",
j_call = "j_call",
junc_len = "junction_length",
nolight = "missing"
)
```

Arguments
-------------------

heavy
:   A tibble containing heavy chain sequences with clone_id

light
:   A tibble containing light chain sequences

nproc
:   number of cores for parallelization

minseq
:   minimum number of sequneces per clone

id
:   name of the column containing sequence identifiers.

seq
:   name of the column containing observed DNA sequences. All 
sequences in this column must be multiple aligned.

clone
:   name of the column containing the identifier for the clone. All 
entries in this column should be identical.

cell_id
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


Details
-------------------

1. Make temporary array containing light chain clones
2. Enumerate all possible V and J combinations
3. Determine which combination is the most frequent
4. Assign sequences with that combination to clone t
5. Copy those sequences to return array
6. Remove all cells with that combination from temp array
7. Repeat 1-5 until temporary array zero.
If there is more than rearrangement with the same V/J
in the same cell, pick the one with the highest non-ambiguous
characters.









