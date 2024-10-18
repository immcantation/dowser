**writeCloneSequences** - *Write the sequences used in tree building to a fasta format. If there are more 
than one tree in airrClone output the sequence id will be followed by "|clone_id".*

Description
--------------------

`writeCloneSequences`Exports the sequences used in tree building.


Usage
--------------------
```
writeCloneSequences(clones, file)
```

Arguments
-------------------

clones
:   tibble `airrClone` objects, the output of 
[formatClones](formatClones.md)

file
:   The file path and name of where the sequences will be saved











