## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
## library(dowser)
## 
## # load example AIRR tsv data
## data(ExampleAirr)
## 
## # Process example data using default settings
## clones = formatClones(ExampleAirr, traits="timepoint", minseq=3)
## 
## # Column shows timepoints in dataset
## print(table(ExampleAirr$timepoint))
## #0   7  14
## #62 102 225
## 
## # Calculate number of tissues sampled in tree
## timepoints = unlist(lapply(clones$data, function(x)
##   length(unique(x@data$timepoint))))
## 
## # Filter to multi-type trees
## clones = clones[timepoints > 1,]
## 
## # Build trees using maximum likelihood (can use alternative builds if desired)
## trees = getTrees(clones, build="pml")
## 


## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
## 
## # correlation test with 10000 repetitions
## test = correlationTest(trees, permutations=10000, time="timepoint")
## print(test)
## 
## # A tibble: 6 × 12
## #  clone_id data       locus  seqs trees      slope     p corre…¹ random…²  min_p
## #     <dbl> <list>     <chr> <int> <list>     <dbl> <dbl>   <dbl>    <dbl>  <dbl>
## #1     3128 <airrClon> N        40 <phylo> -0.00205 0.859  -0.173 -0.0257  0.0667
## #2     3184 <airrClon> N        12 <phylo>  0.00111 0.497   0.649 -0.00429 0.5
## #3     3140 <airrClon> N         9 <phylo>  0.00156 0.335   0.630 -0.00835 0.333
## #4     3192 <airrClon> N         9 <phylo>  0.00739 0.498   0.956 -0.00306 0.5
## #5     3115 <airrClon> N         6 <phylo>  0.00159 0.244   0.565  0.00236 0.25
## #6     3139 <airrClon> N         6 <phylo>  0.00308 0.507   0.821  0.0112  0.5
## 
## 
## # use uniform correlaion test (more sensitive, but higher false positive rate)
## utest = correlationTest(trees, permutations=10000, time="timepoint", perm_type="uniform")
## print(utest)
## 
## # A tibble: 6 × 12
## #  clone_id data       locus  seqs trees      slope      p correlation random_c…¹
## #     <dbl> <list>     <chr> <int> <list>     <dbl>  <dbl>       <dbl>      <dbl>
## #1     3128 <airrClon> N        40 <phylo> -0.00205 0.856       -0.173   0.00146
## #2     3184 <airrClon> N        12 <phylo>  0.00111 0.0768       0.649  -0.00223
## #3     3140 <airrClon> N         9 <phylo>  0.00156 0.114        0.630   0.00205
## #4     3192 <airrClon> N         9 <phylo>  0.00739 0.110        0.956   0.000223
## #5     3115 <airrClon> N         6 <phylo>  0.00159 0.336        0.565   0.00409
## #6     3139 <airrClon> N         6 <phylo>  0.00308 0.172        0.821   0.00431
## 
## 


## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
## library(ggtree)
## 
## # order trees by p value
## test = test[order(test$p),]
## 
## # Plot times on tree with lowest p value (not convincingly evolving)
## plotTrees(test)[[1]] +
##     geom_tippoint(aes(fill=timepoint), pch=21, size=2) +
##     scale_fill_distiller(palette="RdYlBu")


## ---- eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE---------------------
library(dowser)
data(TimeTrees)
TimeTrees = TimeTrees[order(TimeTrees$p),]
plotTrees(TimeTrees)[[1]] + 
    ggtree::geom_tippoint(aes(fill=timepoint), pch=21, size=2) +
    scale_fill_distiller(palette="RdYlBu")

