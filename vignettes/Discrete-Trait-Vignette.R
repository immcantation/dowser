## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
## library(dowser)
## 
## # load example AIRR tsv data
## data(ExampleAirr)
## 
## # trait value of interest
## trait="biopsy"
## 
## # Process example data using default settings
## clones = formatClones(ExampleAirr,
##     traits=trait,num_fields="duplicate_count", minseq=3)
## 
## # Column shows which biopsy the B cell was obtained from
## print(table(ExampleAirr[[trait]]))
## #Lung Nose
## # 145  244
## 
## # Calculate number of tissues sampled in tree
## tissue_types = unlist(lapply(clones$data, function(x)
##   length(unique(x@data[[trait]]))))
## 
## # Filter to multi-type trees
## clones = clones[tissue_types > 1,]
## 
## # Build trees using maximum likelihood (can use alternative builds if desired)
## trees = getTrees(clones, build="pml")
## 


## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
## # the location of the igphyml executable
## # this is location in Docker image, will likely be different if you've set it up yourself
## # note this is the location of the compiled executable file, not just the source folder
## igphyml_location = "/usr/local/share/igphyml/src/igphyml"
## 
## # build trees as before, but use IgPhyML to reconstruct the states of internal
## # nodes using maximum parsimony
## trees = getTrees(clones, build="pml", trait=trait, igphyml=igphyml_location)
## 
## # show internal node (edge) predictions based on maximum parsimony
## plotTrees(trees, tips=trait, nodes=TRUE, node_palette="Set1")[[1]]


## ---- eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE---------------------
library(dowser)
# Load data instead of running phylip
data(BiopsyTrees)
plotTrees(BiopsyTrees, tips="biopsy", nodes=TRUE, node_palette="Set1")[[1]]


## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
## # the location of the igphyml executable
## # this is location in Docker image, will likely be different if you've set it up yourself
## # note this is the location of the compiled executable file, not just the source folder
## igphyml_location = "/usr/local/share/igphyml/src/igphyml"
## 
## # calculate switches along trees compared to 100 random permutations
## # this may take a while, and can be parallelized using nproc
## switches = findSwitches(trees, permutations=100, trait=trait,
##   igphyml=igphyml_location, fixtrees=TRUE)
## 
## # Perform PS test on switches
## ps = testPS(switches$switches)
## print(ps$means)
## # A tibble: 1 x 6
## #  RECON PERMUTE   PLT DELTA STAT   REPS
## #  <dbl>   <dbl> <dbl> <dbl> <chr> <int>
## #1     8     8.6   0.4  -0.6 PS      100
## 
## 
## sp = testSP(switches$switches, alternative="greater")
## print(sp$means)
## # A tibble: 2 x 8
## # Groups:   FROM [2]
## #  FROM  TO    RECON PERMUTE   PGT  DELTA STAT   REPS
## #  <chr> <chr> <dbl>   <dbl> <dbl>  <dbl> <chr> <int>
## #1 Lung  Nose  0.131   0.323     1 -0.192 SP      100
## #2 Nose  Lung  0.869   0.677     0  0.192 SP      100
## 


## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
## # calculate switches along bootstrap distribution of trees
## # build using the 'pml' maximum likelihood option
## # in a real analysis it's important to use at least 100 permutations
## switches = findSwitches(trees, permutations=10, trait=trait,
##   igphyml=igphyml_location, fixtrees=FALSE, build="pml")
## 
## sp = testSP(switches$switches, alternative="greater")
## print(sp$means)
## # A tibble: 2 x 8
## # Groups:   FROM [2]
## #  FROM  TO    RECON PERMUTE   PGT  DELTA STAT   REPS
## #  <chr> <chr> <dbl>   <dbl> <dbl>  <dbl> <chr> <int>
## #1 Lung  Nose  0.168   0.358   1   -0.190 SP       10
## #2 Nose  Lung  0.832   0.642   0.1  0.190 SP       10
## 


## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
## sp = testSP(switches$switches, alternative="greater", permuteAll=TRUE)
## print(sp$means)
## # A tibble: 2 x 8
## # Groups:   FROM [2]
## #  FROM  TO    RECON PERMUTE   PGT   DELTA STAT   REPS
## #  <chr> <chr> <dbl>   <dbl> <dbl>   <dbl> <chr> <int>
## #1 Lung  Nose  0.168   0.241   0.6 -0.0736 SP       10
## #2 Nose  Lung  0.832   0.759   0.4  0.0736 SP       10


## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
## # Downsample each tree to a tip-to-switch ratio of 10 instead of 20
## # this will reduce the false positive rate but also (likely) power
## switches = findSwitches(trees, permutations=100, trait=trait,
##   igphyml=igphyml_location, fixtrees=TRUE, tip_switch=10)
## 
## # didn't have much effect for this dataset
## sp = testSP(switches$switches, alternative="greater")
## print(sp$means)
## # A tibble: 2 x 8
## # Groups:   FROM [2]
## #  FROM  TO    RECON PERMUTE   PGT  DELTA STAT   REPS
## #  <chr> <chr> <dbl>   <dbl> <dbl>  <dbl> <chr> <int>
## #1 Lung  Nose  0.168   0.358   1   -0.190 SP       10
## #2 Nose  Lung  0.832   0.642   0.1  0.190 SP       10
## 


## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
## # the location of the igphyml executable
## # this is location in Docker image, will likely be different if you've set it up yourself
## # note this is the location of the compiled executable file, not just the source folder
## igphyml_location = "/usr/local/share/igphyml/src/igphyml"
## 
## # constant region column name
## trait = "c_call"
## 
## # Vector of human isotypes in the proper order. Isotype switching
## # can only go from left to right (e.g. IGHM to IGHA1). One exception
## # is IGHD, which can switch to IGHM.
## isotypes = c("IGHM","IGHD","IGHG3","IGHG1","IGHA1","IGHG2",
##   "IGHG4","IGHE","IGHA2")
## 
## # Process example data using default settings with "c_call" as a trait value
## clones = formatClones(ExampleAirr, traits=trait, minseq=3)
## 
## # Column shows which constant region associated with a BCR
## print(table(ExampleAirr[[trait]]))
## # IGHA1 IGHA2  IGHD IGHG1 IGHG2 IGHG3 IGHG4  IGHM
## #    55    56    11    58    64    60    63    22
## 
## # Calculate number of istoypes sampled in each tree
## isotype_counts = unlist(lapply(clones$data, function(x)
##   length(unique(x@data[[trait]]))))
## 
## # make model file with irreversibility constraints
## # Will prohibit switches from right to left in the "states" vector
## # IGHD to IGHM switching listed as an exception, since this can occur
## makeModelFile(file="isotype_model.txt", states=isotypes,
##   constraints="irrev", exceptions=c("IGHD,IGHM"))
## 
## # Build trees and predict states at internal nodes using maximum parsimony
## trees = getTrees(clones[isotype_counts > 1,], trait=trait, igphyml=igphyml_location,
##   modelfile="isotype_model.txt", palette="Paired")
## 
## # show internal node (edge) predictions based on maximum parsimony
## plotTrees(trees, tips=trait, nodes=TRUE, node_palette="Paired", ambig="grey")[[1]]


## ---- eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE---------------------
data(IsotypeTrees)
plotTrees(IsotypeTrees, tips="c_call", nodes=TRUE, node_palette="Paired", ambig="grey")[[1]]


## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
## # Downsample each tree to a tip-to-switch ratio of 10 instead of 20
## # this will reduce the false positive rate but also (likely) power
## switches = findSwitches(trees, permutations=100, trait=trait,
##   igphyml=igphyml_location, fixtrees=TRUE, tip_switch=10,
##   modelfile="isotype_model.txt")
## 
## # didn't have much effect for this dataset
## sp = testSP(switches$switches, alternative="greater")
## print(sp$means,n=42)
## # A tibble: 42 × 8
## # Groups:   FROM [7]
## #   FROM  TO      RECON PERMUTE   PGT     DELTA STAT   REPS
## #   <chr> <chr>   <dbl>   <dbl> <dbl>     <dbl> <chr> <int>
## # 1 IGHA1 IGHA2 0.0972  0.0994   0.56 -0.00219  SP      100
## # 2 IGHA1 IGHD  0       0        1     0        SP      100
## # 3 IGHA1 IGHG1 0       0        1     0        SP      100
## # 4 IGHA1 IGHG2 0.00815 0.00833  0.58 -0.000184 SP      100
## # 5 IGHA1 IGHG3 0       0        1     0        SP      100
## # 6 IGHA1 IGHG4 0.00352 0.00409  0.59 -0.000575 SP      100
## # 7 IGHA2 IGHA1 0       0        1     0        SP      100
## # 8 IGHA2 IGHD  0       0        1     0        SP      100
## # 9 IGHA2 IGHG1 0       0        1     0        SP      100
## #10 IGHA2 IGHG2 0       0        1     0        SP      100
## #11 IGHA2 IGHG3 0       0        1     0        SP      100
## #12 IGHA2 IGHG4 0       0        1     0        SP      100
## #13 IGHD  IGHA1 0       0        1     0        SP      100
## #14 IGHD  IGHA2 0       0        1     0        SP      100
## #15 IGHD  IGHG1 0.0144  0.0127   0.42  0.00173  SP      100
## #16 IGHD  IGHG2 0.00926 0.00981  0.54 -0.000545 SP      100
## #17 IGHD  IGHG3 0.0157  0.0113   0.3   0.00440  SP      100
## #18 IGHD  IGHG4 0.00203 0.00271  0.58 -0.000681 SP      100
## #19 IGHG1 IGHA1 0.00632 0.00757  0.59 -0.00126  SP      100
## #20 IGHG1 IGHA2 0.00326 0.00419  0.51 -0.000932 SP      100
## #21 IGHG1 IGHD  0       0        1     0        SP      100
## #22 IGHG1 IGHG2 0.0417  0.0545   0.86 -0.0127   SP      100
## #23 IGHG1 IGHG3 0       0        1     0        SP      100
## #24 IGHG1 IGHG4 0.0451  0.0493   0.66 -0.00427  SP      100
## #25 IGHG2 IGHA1 0       0        1     0        SP      100
## #26 IGHG2 IGHA2 0.00542 0.00511  0.45  0.000316 SP      100
## #27 IGHG2 IGHD  0       0        1     0        SP      100
## #28 IGHG2 IGHG1 0       0        1     0        SP      100
## #29 IGHG2 IGHG3 0       0        1     0        SP      100
## #30 IGHG2 IGHG4 0.0715  0.0581   0.08  0.0133   SP      100
## #31 IGHG3 IGHA1 0.0716  0.0600   0.16  0.0116   SP      100
## #32 IGHG3 IGHA2 0.0508  0.0492   0.37  0.00154  SP      100
## #33 IGHG3 IGHD  0       0        1     0        SP      100
## #34 IGHG3 IGHG1 0.159   0.181    0.96 -0.0225   SP      100
## #35 IGHG3 IGHG2 0.228   0.203    0.05  0.0250   SP      100
## #36 IGHG3 IGHG4 0.161   0.175    0.8  -0.0137   SP      100
## #37 IGHG4 IGHA1 0       0        1     0        SP      100
## #38 IGHG4 IGHA2 0.00596 0.00437  0.27  0.00159  SP      100
## #39 IGHG4 IGHD  0       0        1     0        SP      100
## #40 IGHG4 IGHG1 0       0        1     0        SP      100
## #41 IGHG4 IGHG2 0       0        1     0        SP      100
## #42 IGHG4 IGHG3 0       0        1     0        SP      100


## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
## sp = testSP(switches$switches, alternative="greater", to="IGHA2")
## print(sp$means)
## # A tibble: 8 × 8
## # Groups:   FROM [8]
## #  FROM  TO     RECON PERMUTE   PGT    DELTA STAT   REPS
## #  <chr> <chr>  <dbl>   <dbl> <dbl>    <dbl> <chr> <int>
## #1 IGHA1 IGHA2 0.620   0.634   0.59 -0.0144  SP      100
## #2 IGHD  IGHA2 0       0       1     0       SP      100
## #3 IGHE  IGHA2 0       0       1     0       SP      100
## #4 IGHG1 IGHA2 0.0116  0.0292  0.85 -0.0176  SP      100
## #5 IGHG2 IGHA2 0.0371  0.0252  0.32  0.0119  SP      100
## #6 IGHG3 IGHA2 0.297   0.283   0.45  0.0144  SP      100
## #7 IGHG4 IGHA2 0.0346  0.0290  0.36  0.00563 SP      100
## #8 IGHM  IGHA2 0       0       1     0       SP      100

