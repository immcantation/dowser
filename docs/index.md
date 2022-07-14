Dowser
-------------------------------------------------------------------------------

<img src="img/rod_white_leaf.png" width="200" align=right>

Dowser is part of the [Immcantation](http://immcantation.readthedocs.io) 
analysis framework for Adaptive Immune Receptor Repertoire sequencing 
(AIRR-seq). Dowser provides a set of tools for performing phylogenetic analysis
on B cell receptor repertoires. It supports building and visualizing trees using 
multiple methods, and implements statistical tests for discrete trait analysis
of B cell migration, differentiation, and isotype switching.


Dowser has 6 primary functions:

1. Reconstruct clonal germline sequences.
2. Build B cell lineage trees using multiple methods, such as maximum parsimony, maximum likelihood, and IgPhyML.  
3. Reconstruct intermediate sequences within lineage trees using different methods. 
4. Create publication-quality lineage tree plots.
5. Analyze trees to detect ongoing B cell evolution over time.
6. Understand B cell migration and differentiation (under development).


[What's a dowser?](https://en.wikipedia.org/wiki/Dowsing)

Contact
-------------------------------------------------------------------------------

For help and questions please contact the [Immcantation Group](mailto:immcantation@googlegroups.com)
or use the [issue tracker](https://bitbucket.org/kleinstein/dowser/issues?status=new&status=open).



# Dependencies

**Depends:** ggplot2  
**Imports:** alakazam, ape, Biostrings, dplyr, ggtree, graphics, gridExtra, markdown, methods, phangorn, phylotate, RColorBrewer, rlang, shazam, stats, stringr, tidyselect, tidyr, utils  
**Suggests:** knitr, rmarkdown, testthat


# Authors

[Kenneth Hoehn](mailto:kenneth.hoehn@yale.edu) (aut, cre)  
[Susanna Marquez](mailto:susanna.marquez@yale.edu) (ctb)  
[Jason Vander Heiden](mailto:jason.vanderheiden@gmail.com) (ctb)  
[Steven Kleinstein](mailto:steven.kleinstein@yale.edu) (aut, cph)


# Citing


To cite the dowser package in publications, please use:

Hoehn K, Pybus O, Kleinstein S (2020). “Phylogenetic analysis of
migration, differentiation, and class switching in B cells.” _bioRxiv_.
doi: 10.1101/2020.05.30.124446 (URL:
https://doi.org/10.1101/2020.05.30.124446), <URL:
https://doi.org/10.1101/2020.05.30.124446>.

A BibTeX entry for LaTeX users is

  @Article{,
    style = {citation},
    title = {Phylogenetic analysis of migration, differentiation, and class switching in B cells.},
    author = {Kenneth B. Hoehn and Oliver G. Pybus and Steven H. Kleinstein},
    year = {2020},
    journal = {bioRxiv},
    url = {https://doi.org/10.1101/2020.05.30.124446},
    doi = {10.1101/2020.05.30.124446},
  }



# License

AGPL-3
