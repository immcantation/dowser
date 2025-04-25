Version 2.3.2 April 25, 2025
-------------------------------------------------------------------------------

+ Added light chain UCA inference to getTreesAndUCAs and addressed minor bugs

Version 2.3.1 March 13, 2025
-------------------------------------------------------------------------------

+ Added a UCA inference function called getTreesAndUCAs

Version 2.3 October 18, 2024
-------------------------------------------------------------------------------

+ Fixed bugs in and added features to buildRAxML
+ Fixed bugs in resolveLightChains
+ Added trimming option to createGermlines

Version 2.2.1 May 23, 2024
-------------------------------------------------------------------------------

+ Updated correlationTest
+ Fixed bugs in buildRAxML
+ Added ASR check for version changes in phangorn

Version 2.2.0 May 7, 2024
-------------------------------------------------------------------------------

+ added suggestion for pwalign dependency
+ set filter_stop=FALSE in formatClones for better speed

Version 2.1.0 December 21, 2023
-------------------------------------------------------------------------------

+ plotTrees uses palette instead of tip_palette/node_palette
+ palette options deprecated in getTrees and sub-functions
+ getPalette colors Germline as black by default
+ ambiguous nodes are now grey by default

Version 2.0.1: October 25, 2023
-------------------------------------------------------------------------------

+ Fixed bugs in buildRAxML and getTrees
+ Added an option to formatClones to consider light chain traits

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
