**checkDivergence** - *Compare divergence along a tree in terms of mutations (sum of branches)
for each tip and reconstructed internal node 
to its Hamming distance from the germline. Divergence should never be less than Hamming distance. 
A threshold of -1 is used to represent 
1 full mutation difference. The function will throw a warning if any trees
cross this threshold*

Description
--------------------

Compare divergence along a tree in terms of mutations (sum of branches)
for each tip and reconstructed internal node 
to its Hamming distance from the germline. Divergence should never be less than Hamming distance. 
A threshold of -1 is used to represent 
1 full mutation difference. The function will throw a warning if any trees
cross this threshold


Usage
--------------------
```
checkDivergence(clones, threshold = -1, verbose = TRUE, germline = "Germline")
```

Arguments
-------------------

clones
:   a tibble of clones and trees, output from [getTrees](getTrees.md)

threshold
:   Minimum allowed value of divergence minus Hamming distance

verbose
:   Print whether all trees passed

germline
:   ID of the tree's predicted germline sequence




Value
-------------------

`tibble` showing the clone_id, sequence_id, as well as tree-based
divergence, hamming distance, and difference between the two.









