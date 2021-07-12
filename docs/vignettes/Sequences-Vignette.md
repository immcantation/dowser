# Reconstructing intermediate sequences

Dowser automatically reconstructs intermediate sequences as part of the `getTrees` function. These are stored in the `nodes` list contained in each `phylo` object.

## Get intermediate sequences

First, collapse internal nodes with identical sequences. This will significantly clean up the visualization. Then, visualize the trees using `plotTrees` but with the node_nums parameter set. This will display the number of each internal node.

To obtain the IMGT-gapped sequence for each reconstructed node, specify the clone ID and node number in the `getSeq` function.

```r

# Collapse nodes with identical sequences. This will 
clones = collapseNodes(clones)

plots = plotTrees(clones, tips="c_call", tipsize=2,
    node_nums=TRUE, labelsize=7)

treesToPDF(plots[[1]], file="trees.pdf", nrow=1, ncol=1)

getSeq(clones, node=54, clone=3128)

# Output:
#                                                                      N 
# "GAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTACAGCCAGGGCGGTCCCTGAGACTCTCCTGTACAGCTTCTGGATTCACCTTT............AGTGAYTATGCTATGAGCTGGTTCCGCCAGGCTCCAGGGAAGGGKCTGGAGTGGGTGGGTTTCATTAGAAGCAGACGTTTTGGTGGGACGCCGGACTACGCCGCGTCAGTGAGA...GACAGATTCACCATTTCAAGAGACGATTCCAAAAGCATCGCCTATCTGCAAATGAACAGCCTGAAAACCGAGGACACAGCCGTGTATTTTTGTAGTAGAGATCTCGCGGTTATATCCACAATAGCTGGTACTAACTGGTTCGACCCCAGGGGCCAGGGAGCCCTGGTCACCGTCTCCTCAGNN" 
```
![Tree with nodes labeled](figures/Sequences-Vignette-all.png)