**makeModelFile** - *Make a parsimony model file*

Description
--------------------

`makeModelFile` Filler


Usage
--------------------
```
makeModelFile(file, states, constraints = NULL, exceptions = NULL)
```

Arguments
-------------------

file
:   model file name to write.

states
:   vector of states to include in model.

constraints
:   constraints to add to model.

exceptions
:   vector of comma-separated states that are 
exceptions to constraints




Value
-------------------

Name of model file


Details
-------------------

Currently the only option for `constraints` is "irrev", which
forbids switches moving from left to right in the `states` vector.




See also
-------------------

[readModelFile](readModelFile.md), [getTrees](getTrees.md), [findSwitches](findSwitches.md)






