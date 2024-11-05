**getAllSeqs** - *Return all tip and internal node sequences*

Description
--------------------

`getNodeSeq` Sequence retrieval function.


Usage
--------------------
```
getAllSeqs(data, imgt_gaps = TRUE)
```

Arguments
-------------------

data
:   a tibble of `airrClone` objects with reconstructed trees, 
the output of [getTrees](getTrees.md)

imgt_gaps
:   include a column of gapped sequences?




Value
-------------------

A tibble with sequence information for each tip and internal node
of a set of trees.


Details
-------------------

Column names:
clone_id = clone id 
node_id = name of node, either the sequence name if a tip or Node<number> if internal node
node = node number in tree. Tips are nodes 1:<number of tips>.
locus = locus of sequence
sequence = ungapped sequence, either observed for tips or reconstructed for internal nodes
sequence_alignment = sequence with IMGT gaps (optional)




See also
-------------------

[getTrees](getTrees.md) [getNodeSeq](getNodeSeq.md)






