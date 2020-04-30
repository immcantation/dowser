**ExampleTrees** - *Example Ig lineage trees*

Description
--------------------

A set of Ig lineage trees generated from the `ExampleDb` file, subset to
only those trees with at least four nodes.


Usage
--------------------
```
ExampleTrees
```




Format
-------------------

A list of igraph objects output by buildPhylipLineage.
Each node of each tree has the following annotations (vertex attributes):

+ `sample`:    Sample identifier(s). Time in relation to vaccination.
+ `isotype`:   Isotype assignment(s). 
+ `duplication_count`:  Copy count (number of duplicates) of the sequence.





See also
-------------------

[ExampleTrees](ExampleTrees.md)






