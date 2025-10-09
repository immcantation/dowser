**plotSkylines** - *Simple function for plotting Bayesian skyline plots*

Description
--------------------

`plotSkylines` Simple Bayesian skyline plots


Usage
--------------------
```
plotSkylines(clones, file = NULL, width = 8.5, height = 11)
```

Arguments
-------------------

clones
:   output from getTrees using BEAST

file
:   pdf file name for printing plots

width
:   width of plot in inches if file specified

height
:   height of plot in inches if file specified

...
:   optional arguments passed to grDevices::pdf




Value
-------------------

if no file specified, a list of ggplot objects. If file specified
will plot to specified file




See also
-------------------

[getSkylines](getSkylines.md) [readBEAST](readBEAST.md) [getTrees](getTrees.md)






