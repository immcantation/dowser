**create_starting_tree** - *Takes an airr clone object and tree and returns BEAST2 XML for setting the starting tree*

Description
--------------------

Takes an airr clone object and tree and returns BEAST2 XML for setting the starting tree


Usage
--------------------
```
create_starting_tree(
clone,
id,
tree,
include_germline_as_tip,
tree_states,
start_edge_length
)
```

Arguments
-------------------

clone
:   an `airrClone` object

id
:   unique identifer for this analysis

tree
:   starting tree, either a phylo object or a newick string

include_germline_as_tip
:   include the germline as a tip

tree_states
:   use states in the starting tree?

start_edge_length
:   edge length to use for all branches in starting tree




Value
-------------------

String of XML setting the starting tree









