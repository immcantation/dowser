Dowser
-------------------------------------------------------------------------------

Dowser is part of the [Immcantation](http://immcantation.readthedocs.io) 
analysis framework for Adaptive Immune Receptor Repertoire sequencing 
(AIRR-seq) and provides a set of tools for performing phylogenetic analysis
on B cell receptor repertoires

Notice
-------------------------------------------------------------------------------

Dowser isn’t officially released yet, so please let us know if you’re planning to publish anything using it. A lot of the features are still in active development and will likely change in the near future. Feedback is be greatly appreciated.

Install
-------------------------------------------------------------------------------

Currently dowser can only be built from the [source code](http://bitbucket.org/kleinstein/dowser),
first install the build dependencies:

```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "Rcpp"))
```

To install the latest development code via devtools, along with a development branch of alakazam:

```R
library(devtools)
install_bitbucket("kleinstein/alakazam@dowser")
install_bitbucket("kleinstein/dowser@master")
```

Quick start: lineage tree reconstruction and discrete trait analysis
----------------------------------------------------------------------------

The following commands, entered into an R terminal, go through basic operations of building B cell lineage trees and discrete trait analysis using the PS, SC, and SP tests. By default, trees are built using maximum parsimony (using the `pratchet` function of the R package `phangorn`). Switches among trait values are calculated using IgPhyML (>v1.1.0). You must supply the path to the compiled IgPhyML executible for function `bootstrapTrees` to work. See https://igphyml.readthedocs.io for IgPhyML installation details.

```R
# Load required packages
library(alakazam)
library(dowser)

# load example AIRR tsv data
data(ExampleDb)

# Process example data into proper format
clones <- formatClones(ExampleDb, trait="c_call")

# Build maxmimum parsimony trees for first two clones using 
# phangorn package in R
trees <- getTrees(clones[1:2,])

# simple tree plotting with ggtree R package
plots <- plotTrees(trees, tips="c_call")

# plot tree of first clone
plots[[1]]

# path to igphyml executible
igphyml <- "PATH/TO/IGPHYML"

# Count switches across bootstrap distribution of trees
bootstraps <- bootstrapTrees(clones[1:2,], trait="c_call",
	bootstraps=10, igphyml=igphyml)

# perform PS, SC, and SP tests on switch counts
testPS(bootstraps$switches)
testSC(bootstraps$switches)
testSP(bootstraps$switches)

```
