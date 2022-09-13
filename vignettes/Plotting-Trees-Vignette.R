## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------

library(dowser)

data(ExampleClones)

ExampleClones = ExampleClones[1:2,]

plots = plotTrees(ExampleClones)

#Plot the largest tree
#To plot second largest tree, use plots[[2]], and so on
plots[[1]]



## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------

# Re-scale branches to represent mutations between nodes
ExampleClones_m = scaleBranches(ExampleClones, edge_type="mutations")

# Plot, set scale bar to represent 10 mutations
plots = plotTrees(ExampleClones_m, scale=10)

#Plot the largest tree
plots[[1]]



## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Plot tree with sequence isotype at the tips.
plots = plotTrees(ExampleClones, tips="c_call")

#Plot the largest tree
plots[[1]]


## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Plot tree with sequence isotype at the tips, with sizes set to number of duplicates
plots = plotTrees(ExampleClones, tips="c_call", tipsize="duplicate_count")

#Plot the largest tree
plots[[1]]


## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# These calls create the same plot:

# Plot tree with sequence isotype at the tips, with palette "Set1"
plots = plotTrees(ExampleClones, tips="c_call", tipsize=2,
    tip_palette="Paired")

# or, specify a named palette vector
custom_palette=c(
    "IGHM"    ="#a6cee3",
    "IGHD"    ="#1f78b4",
    "IGHG3"   ="#b2df8a",
    "IGHG1"   ="#33a02c",
    "IGHA1"   ="#fb9a99",
    "IGHG2"   ="#e31a1c",
    "IGHG4"   ="#fdbf6f",
    "IGHE"    ="#ff7f00",
    "IGHA2"   ="#cab2d6",
    "Germline"="#6a3d9a")

plots = plotTrees(ExampleClones, tips="c_call", tipsize=2,
    tip_palette=custom_palette)

# or, use the getPalette function to create a named palette vector
custom_palette = getPalette(c("IGHM","IGHD","IGHG3","IGHG1","IGHA1",
    "IGHG2","IGHG4","IGHE","IGHA2","Germline"), "Paired")

plots = plotTrees(ExampleClones, tips="c_call", tipsize=2,
    tip_palette=custom_palette)

#Plot the largest tree
plots[[1]]


## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
library(ggtree)

plots = plotTrees(ExampleClones, tips="c_call", tipsize=2)

#Plot the largest tree
treeplot = plots[[1]] + geom_tiplab() + 
    geom_vline(xintercept=c(0.05,0.1,0.15,0.2,0.25),
        linetype="dashed",color="grey") + xlim(0,0.3) +
    ggtitle("Example B cell tree")

treeplot


## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
## plots = plotTrees(ExampleClones, tips="c_call", tipsize=2)
## 
## # you can also pass arguments you would pass to grDevices::pdf, like width and height
## # here, we plot 4 trees per page (2 rows, 2 columns)
## treesToPDF(plots, file="trees.pdf", nrow=2, ncol=2)

