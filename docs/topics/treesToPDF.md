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
:   size of tip shape objects

...
:   optional arguments passed to grDevices::pdf




Value
-------------------

a PDF of tree plots



Examples
-------------------

```R
### Not run:
data(ExampleDb)
# ExampleDb$sample_id <- sample(ExampleDb$sample_id)
# clones <- formatClones(ExampleDb, trait="sample_id")
# trees <- getTrees(clones)
# plots <- plotTrees(trees)
# treesToPDF(plots,"test.pdf",width=5,height=6)
```



See also
-------------------

[plotTrees](plotTrees.md)






