# Quickstart

Quick start of lineage tree reconstruction.

```r
# Load required packages
library(alakazam)
library(dowser)

# load example AIRR tsv data
data(ExampleDb)

# Process example data into proper format
clones = formatClones(ExampleDb)

# Build maxmimum parsimony trees for first two clones using 
# phangorn package in R
trees = getTrees(clones[1:2,])

# simple tree plotting with ggtree R package
plots = plotTrees(trees)

# plot tree of first clone
plots[[1]]
```
