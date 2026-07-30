**getSeqPath** - *Return all sequences along the tree from the germline to a specified node*

Description
--------------------

Return all sequences along the tree from the germline to a specified node


Usage
--------------------
```
getSeqPath(
data,
node,
tree = NULL,
clone = NULL,
gaps = TRUE,
translate = FALSE,
rm_x_aa = FALSE,
dot_notation = FALSE,
verbose = FALSE
)
```

Arguments
-------------------

data
:   a tibble of `airrClone` objects, the output of 
[getTrees](getTrees.md)

node
:   numeric node in tree (see details)

tree
:   a `phylo` tree object containing `node`

clone
:   if `tree` not specified, supply clone ID in `data`

gaps
:   add IMGT gaps to output sequences?

translate
:   Tranlate into amino acids?

rm_x_aa
:   Remove amino acids translated to X

dot_notation
:   Represent sites with germline character as "." and ambiguous sites as "?"

verbose
:   print out extra information?




Value
-------------------

A data table where each row shows a node/node id, with that node's sequence(s)
in the remaining columns with respective locus labels.




See also
-------------------

[getNodeSeq](getNodeSeq.md) [getAllSeqs](getAllSeqs.md) [dfToFasta](dfToFasta.md)






