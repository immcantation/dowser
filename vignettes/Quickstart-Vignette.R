## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Load required packages
library(alakazam)
library(dowser)

# load example AIRR tsv data
data(ExampleAirr)

# subset data for this example
ExampleAirr = ExampleAirr[ExampleAirr$clone_id %in% c("3170", "3184"),]

# Process example data into proper format, store isotype (optional)
clones = formatClones(ExampleAirr, trait="c_call")

# Build maxmimum parsimony trees for first two clones using 
# phangorn package in R
trees <- getTrees(clones)

# simple tree plotting with ggtree R package with isotypes at tips
plots <- plotTrees(trees, tips="c_call",tipsize=2)

# plot tree of largest clone
plots[[1]]

