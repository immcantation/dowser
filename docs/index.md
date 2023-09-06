Dowser
-------------------------------------------------------------------------------

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
6. Understand B cell migration and differentiation.


[What's a dowser?](https://en.wikipedia.org/wiki/Dowsing)

Contact
-------------------------------------------------------------------------------

For help and questions please contact the [Immcantation Group](mailto:immcantation@googlegroups.com)
or use the [issue tracker](https://bitbucket.org/kleinstein/dowser/issues?status=new&status=open).



## Dependencies

**Depends:** ggplot2  
**Imports:** alakazam, ape, Biostrings, dplyr, ggtree, graphics, gridExtra, markdown, methods, phangorn, phylotate, RColorBrewer, rlang, shazam, stats, stringr, tidyselect, tidyr, utils  
**Suggests:** knitr, rmarkdown, testthat


## Authors

[Kenneth Hoehn](mailto:kenneth.hoehn@yale.edu) (aut, cre)  
[Cole Jensen](mailto:cole.jensen@yale.edu) (ctb)  
[Susanna Marquez](mailto:susanna.marquez@yale.edu) (ctb)  
[Jason Vander Heiden](mailto:jason.vanderheiden@gmail.com) (ctb)  
[Steven Kleinstein](mailto:steven.kleinstein@yale.edu) (aut, cph)


## Citing

To cite the dowser package in publications, please use

  Hoehn K, Pybus O, Kleinstein S (2022). “Phylogenetic analysis of
  migration, differentiation, and class switching in B cells.” _PLoS
  Computational Biology_. doi:10.1371/journal.pcbi.1009885
  <https://doi.org/10.1371/journal.pcbi.1009885>,
  <https://doi.org/10.1371/journal.pcbi.1009885>.

If you use the [correlationTest](vignettes/Measurable-Evolution.md)
function for measurable evolution, please also cite

  Hoehn K, Turner J, Miller F, Jiang R, Ellebedy A, Pybus O, Kleinstein
  S (2021). “Human B cell lineages associated with germinal centers
  following influenza vaccination are measurably evolving.” _eLife_.
  doi:10.7554/eLife.70873 <https://doi.org/10.7554/eLife.70873>,
  <https://elifesciences.org/articles/70873>.

If you additionally use [IgPhyML](https://igphyml.readthedocs.io) for
building trees please also cite

  Hoehn K, Van der Heiden J, Zhou J, Lunter G, Pybus O, Kleinstein S
  (2019). “Repertoire-wide phylogenetic models of B cell molecular
  evolution reveal evolutionary signatures of aging and vaccination.”
  _PNAS_. doi:10.1073/pnas.1906020116
  <https://doi.org/10.1073/pnas.1906020116>,
  <https://doi.org/10.1073/pnas.1906020116>.



## License

AGPL-3
