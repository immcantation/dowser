**treesToPDF** - *Simple function for plotting a lot of trees into a pdf*

Description
--------------------

`treesToPDF` exports trees to a pdf in an orderly fashion


Usage
--------------------
```
treesToPDF(plots, file, nrow = 2, ncol = 2, ...)
```

Arguments
-------------------

plots
:   list of tree plots (from plotTrees)

file
:   output file name

nrow
:   number of rows per page

ncol
:   number of columns per page

...
:   optional arguments passed to grDevices::pdf




Value
-------------------

a PDF of tree plots



Examples
-------------------

```R
### Not run:
data(ExampleClones)
# trees <- getTrees(ExampleClones[10,])
# plots <- plotTrees(trees)
# treesToPDF(plots,"test.pdf",width=5,height=6)

```



See also
-------------------

[plotTrees](plotTrees.md)






