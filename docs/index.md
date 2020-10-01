Dowser
-------------------------------------------------------------------------------

Dowser is part of the [Immcantation](http://immcantation.readthedocs.io) 
analysis framework for Adaptive Immune Receptor Repertoire sequencing 
(AIRR-seq). Dowser provides a set of tools for performing phylogenetic analysis
on B cell receptor repertoires. It supports building and visualizing trees using 
multiple methods, and implements statistical tests for discrete trait analysis
of B cell migration, differentiation, and isotype switching.

Dowser is released under the AGPL-3 license.


Notice
-------------------------------------------------------------------------------

Dowser isn’t officially released yet, so please let us know if you’re planning to publish anything using it. All features are still in active development and will likely change in the near future. Feedback is be greatly appreciated.

Install
-------------------------------------------------------------------------------

Dowser can be built from the [source code](http://bitbucket.org/kleinstein/dowser), or used in the Immcantation Docker container development build. The Docker container currently the only way to use the igphyml dependent features of Dowser on Windows.

To use Dowser in the Docker image, simply pull and run the immcantation/suite:devel [Docker image](https://immcantation.readthedocs.io/en/stable/docker/intro.html). Instructions in hyperlink.

Alternatively, to install the source code, first install the build dependencies:

```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "Rcpp"))
```

To install the latest development code via devtools, along with the development version of Alakazam:

```R
library(devtools)
install_bitbucket("kleinstein/alakazam@master")
install_bitbucket("kleinstein/dowser@master")
```

Quick start: lineage tree reconstruction and discrete trait analysis
----------------------------------------------------------------------------

The following commands, entered into an R terminal, go through basic operations of building B cell lineage trees and discrete trait analysis using the PS, SC, and SP tests. By default, trees are built using maximum parsimony (using the `pratchet` function of the R package `phangorn`). Switches among trait values are calculated using [IgPhyML](https://igphyml.readthedocs.io) (>v1.1.0). You must supply the path to the compiled IgPhyML executible for function `bootstrapTrees` to work. If using the Docker container, IgPhyML will already be installed. If using the source code version, see https://igphyml.readthedocs.io for IgPhyML installation details.

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

# Path to igphyml executible
# The path in Docker image is supplied. If installed from 
# source this will likely be different.
igphyml <- "/usr/local/share/igphyml/src/igphyml"

# Count switches across bootstrap distribution of trees
bootstraps <- bootstrapTrees(clones[1:2,], trait="c_call",
	bootstraps=10, igphyml=igphyml)

# perform PS, SC, and SP tests on switch counts
testPS(bootstraps$switches)
testSC(bootstraps$switches)
testSP(bootstraps$switches)

```


# Dependencies

**Depends:** ggplot2  
**Imports:** alakazam, ape, dplyr, ggtree, graphics, grid, gridExtra, igraph, lazyeval, Matrix, methods, phangorn, phylotate, progress, RColorBrewer, Rcpp, readr, rlang, scales, seqinr, stats, stringi, stringr, tibble, tidyselect, tidyr, utils  
**Suggests:** knitr, rmarkdown, testthat  
**Extends:** FALSE


# Authors

[Kenneth Hoehn](mailto:kenneth.hoehn@yale.edu) (aut, cre)  
[Steven Kleinstein](mailto:steven.kleinstein@yale.edu) (aut, cph)


# Citing


To cite the dowser package in publications, please use:

Hoehn K, Pybus O, Kleinstein S (2020). “Phylogenetic analysis of
migration, differentiation, and class switching in B cells.”
_bioRxiv_, doi: https://doi.org/10.1101/2020.05.30.124446.

A BibTeX entry for LaTeX users is

  @Article{,
    style = {citation},
    title = {Phylogenetic analysis of migration, differentiation, and class switching in B cells.},
    author = {Kenneth B. Hoehn and Oliver G. Pybus and Steven H. Kleinstein},
    year = {2020},
    journal = {bioRxiv},
    DOI = {https://doi.org/10.1101/2020.05.30.124446},
  }

