# Release Notes

Version 2.0.0: September 8, 2023
-------------------------------------------------------------------------------

+ Added paired heavy and light chain single cell integration into formatClones
+ Added lone light chain integration into formatClones
+ Added resolveLightChains, a function that assigned clones 'clonal_subgroups'
+ Added a vignette about how to create H+L trees

Version 1.2.0: May 30, 2023
-------------------------------------------------------------------------------

+ Fixed bugs in createGermlines, plotTrees, and getTrees
+ Added getBootstraps, a bootstrapping function separate from findSwitches
+ Added RAxML build option
+ Added partitioned IgPhyML and RAxML build options (experimental)
+ Added calcRF for RF distance calculations.


Version 1.1.0:  October 4, 2022
-------------------------------------------------------------------------------

+ createGermlines now works on all loci simultaneously
+ Vignettes for measurable evolution and discrete trait analysis
+ Corrected bugs in createGermlines and findSwitches

Version 1.0.0:  April 8, 2022
-------------------------------------------------------------------------------

+ Deprecated bootstrapTrees in favor of findSwitches
+ Deprecated getSeq in favor of getNodeSeq
+ Added vignettes for discrete trait analysis and measurable evolution
+ Changed permutation to perm_type in correlationTest
+ Added common_scale to plotTrees
+ Removed trace and warning printouts from getTrees(build="pml")


Version 0.1.0:  July 13, 2021
-------------------------------------------------------------------------------

+ Initial CRAN release.


Version 0.0.3:  December 12, 2020
-------------------------------------------------------------------------------

+ Correlation test added.
+ Additional vignettes added.


Version 0.0.2:  August 28, 2020
-------------------------------------------------------------------------------

+ Integration with new Alakazam lineage functions


Version 0.0.0.999:  April 14, 2020
-------------------------------------------------------------------------------

+ Base functions getTrees, bootstrapTrees, and plotTrees
