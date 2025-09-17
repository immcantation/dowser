## ----prepare-data, eval=FALSE------------------------------------------------------------------------------------------------
## library(dowser)
## library(dplyr)
## library(ggtree)
## 
## # load example AIRR tsv data
## data(ExampleAirrTyCHE)
## 
## # set up time/date trait
## ExampleAirrTyche$sample_time <- as.numeric(ExampleAirrTyche$sample_time)
## 
## # trait value of interest
## trait="location"
## 
## clones <- formatClones(
##   ExampleAirrTyche,
##   traits = c(trait, "sample_time"),
##   germ   = "germline_alignment"
## )
## 
## # Column shows which location the B cell was obtained from
## print(table(ExampleAirrTyCHE[[trait]]))


## ----print-table, echo=FALSE, warning=FALSE, message=FALSE-------------------------------------------------------------------
library(dowser)
data(ExampleAirrTyCHE)
trait="location"

print(table(ExampleAirrTyCHE[[trait]]))


## ----estimate-gc-clock-rate, eval=FALSE--------------------------------------------------------------------------------------
## gc_cells = filter(ExampleAirrTyCHE, location=="germinal_center")
## gcf = formatClones(gc_cells, traits=c("location","sample_time"),
## 	germ="germline_alignment")
## 
## 
## gctrees = getTrees(gcf, build="pml", sub_model="HKY")
## 
## plotTrees(gctrees)[[1]] + geom_tippoint(aes(color=sample_time))


## ----eval=FALSE--------------------------------------------------------------------------------------------------------------
## # edit to your BEAST installation path
## beast <- "/Applications/BEAST 2.7.7/bin/"
## 
## # estimate clock rate of GC B cells
## # Note, this will take several minutes
## # if you don't care about convergence, reduce mcmc_length
## gctree = getTimeTreesIterate(gcf,
## 	beast=beast,
## 	template="StrictClock_Standard_EmpFreq.xml",
## 	dir="temp",
## 	id="gc_strict",
## 	time="sample_time",
## 	mcmc_length=1e7,
## 	iterations=5,
## 	nproc=5,
## 	CLOCK_RATE_INIT=0.001,
## 	KAPPA_PRIOR_M=0.67,
## 	KAPPA_PRIOR_S=0.2)
## 
## 
## gcrate_tree = mean(sapply(gctree$parameters, function(x)filter(x,item=="geneticClockRate")$mean))
## print(gcrate_tree)

## ----echo=FALSE--------------------------------------------------------------------------------------------------------------
print(0.0004308)


## ----eval=FALSE--------------------------------------------------------------------------------------------------------------
## gcrate_slope = mean(correlationTest(gctrees, time="sample_time")$slope)
## print(gcrate_slope)


## ----echo=FALSE--------------------------------------------------------------------------------------------------------------
print(0.0003716162)


## ----eval=FALSE--------------------------------------------------------------------------------------------------------------
## mixed_trees <- getTimeTreesIterate(
##   clones,
##   beast    = beast,
##   template = "TraitLinkedExpectedOccupancy_EstTraitClockRates_EmpFreq.xml",
##   trait    = trait,
##   time     = "sample_time",
##   dir      = "temp",
##   id       = "tyche_eo_est",
##   log_every = "auto",
##   nproc     = 5,
##   KAPPA_PRIOR_M = 0.67,
##   KAPPA_PRIOR_S = 0.2,
##   TRAIT_RATE_MEAN_1 = gcrate_tree,
##   TRAIT_RATE_MEAN_2 = 0.000001,
##   TRAIT_RATE_SIGMA_1 = gcrate_tree * 0.01,
##   TRAIT_RATE_SIGMA_2 = 0.00001,
##   RATE_INDICATORS = "1 0",
##   TRANSITION_RATE_ALPHA_1 = 0.1,
##   TRANSITION_RATE_ALPHA_2 = 1.0,
##   TRANSITION_RATE_BETA_1  = 0.1,
##   TRANSITION_RATE_BETA_2  = 1.0,
##   log_target   = 2000,
##   mcmc_length  = 1e7,
##   ignore       = c("freqParameter"),
##   iterations   = 10
## )


## ----plot-mixed-trees-result, eval=FALSE-------------------------------------------------------------------------------------
## plotTrees(mixed_trees)[[1]] + geom_point(aes(color=location))


## ----read-beast-output, eval=FALSE-------------------------------------------------------------------------------------------
## mixed_trees <- readBEAST(clones, dir="temp", id="tyche_eo_est", beast=beast, trait=trait)

