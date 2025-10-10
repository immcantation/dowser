# [![](http://cranlogs.r-pkg.org/badges/grand-total/dowser)](https://www.r-pkg.org/pkg/dowser) [![](https://cranlogs.r-pkg.org/badges/dowser)](https://www.r-pkg.org/pkg/dowser) [![](https://img.shields.io/static/v1?label=AIRR-C%20sw-tools%20v1&message=compliant&color=008AFF&labelColor=000000&style=plastic)](https://docs.airr-community.org/en/stable/swtools/airr_swtools_standard.html)

**IMPORTANT!** 
Dowser's source code has moved to https://github.com/immcantation/dowser

To update Git configuration settings use:

```
   git config user.email "your-gh-user@email.com"
   git config user.name "your-gh-user-name"
   git remote set-url origin git@github.com:immcantation/dowser.git
```

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


Documentation and tutorials
-------------------------------------------------------------------------------

Full documentation, including tutorials and vignettes: https://dowser.readthedocs.io


Contact
-------------------------------------------------------------------------------

If you need help or have any questions, please contact the [Immcantation Group](mailto:immcantation@googlegroups.com).

If you have discovered a bug or have a feature request, you can open an issue using the [issue tracker](https://github.com/immcantation/dowser/issues).

To receive alerts about Immcantation releases, news, events, and tutorials, join the [Immcantation News](https://groups.google.com/g/immcantation-news) Google Group. [Membership settings](https://groups.google.com/g/immcantation-news/membership) can be adjusted to change the frequency of email updates.



## Dependencies

**Depends:** ggplot2  
**Imports:** airr, alakazam, ape, Biostrings, dplyr, ggtree, graphics, gridExtra, markdown, methods, phangorn, phylotate, RColorBrewer, rlang, shazam, stats, stringr, tidyselect, tidyr, utils, treeio  
**Suggests:** knitr, rmarkdown, testthat, pwalign, BiocManager


## Authors

[Kenneth Hoehn](mailto:kenneth.b.hoehn@dartmouth.edu) (aut, cre)  
[Cole Jensen](mailto:cole.jensen@yale.edu) (aut)  
[Jessie Fielding](mailto:jessie.jo.fielding@dartmouth.edu) (aut)  
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

If you construct [paired heavy and light chain
trees](vignettes/Resolve-Light-Chains-Vignette.md), please also cite

  Jensen C, Sumner J, Kleinstein S, Hoehn K (2024). “Inferring B Cell
  Phylogenies from Paired H and L Chain BCR Sequences with Dowser.”
  _The Journal of Immunology_. doi:10.4049/jimmunol.2300851
  <https://doi.org/10.4049/jimmunol.2300851>,
  <https://doi.org/10.4049/jimmunol.2300851>.

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
