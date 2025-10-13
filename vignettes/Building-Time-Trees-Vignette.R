## # Choose appropriate version for your architecture (x86 or aarch64)
## BEAST=BEAST.v2.7.7.Linux.x86.tgz # or BEAST=BEAST.v2.7.7.Linux.aarch64.tgz
## 
## # download file and uncompress
## curl -O https://github.com/CompEvol/beast2/releases/download/v2.7.7/$BEAST
## tar -xvzf $BEAST
## 
## # optionally remove the compressed file
## rm $BEAST
## 
## # run BEAST, at least with help, to allow it to set up its directories
## ~/beast/bin/beast -help
## 
## # install BEAST Classic package
## ~/beast/bin/packagemanager -add BEAST_CLASSIC
## 
## # install TyCHE package (currently not released on BEAST package manager)
## curl -O https://github.com/hoehnlab/tyche/releases/download/v0.0.3/TyCHE.v0.0.3.zip
## unzip -o -d ~/.beast/2.7/TyCHE TyCHE.v0.0.3.zip
## rm -f TyCHE.v0.0.3.zip
## 
## # install TyCHE package once it is released on BEAST package manager
## # ~/beast/bin/packagemanager -add TyCHE
## 
## # install rootfreqs package
## ROOTFREQS=rootfreqs.package.v0.0.2.zip
## curl -O https://github.com/rbouckaert/rootfreqs/releases/download/v0.0.2/$ROOTFREQS
## unzip -o -d ~/.beast/2.7/rootfreqs $ROOTFREQS
## rm -f $ROOTFREQS
## 

## ----prepare-data, eval=FALSE-------------------------------------------------------------------------------------------------
## library(dowser)
## library(dplyr)
## library(ggtree)
## 
## # load example AIRR tsv data
## data(ExampleAirrTyCHE)
## 
## # set up time/date trait
## ExampleAirrTyCHE$sample_time <- as.numeric(ExampleAirrTyCHE$sample_time)
## 
## # trait value of interest
## trait="location"
## 
## clones <- formatClones(
##   ExampleAirrTyCHE,
##   traits = c(trait, "sample_time"),
##   germ   = "germline_alignment"
## )
## 
## # Column shows which location the B cell was obtained from
## print(table(ExampleAirrTyCHE[[trait]]))


## ----print-table, echo=FALSE, warning=FALSE, message=FALSE--------------------------------------------------------------------
library(dowser)
data(ExampleAirrTyCHE)
trait="location"

print(table(ExampleAirrTyCHE[[trait]]))


## ----estimate-gc-clock-rate, eval=FALSE---------------------------------------------------------------------------------------
## gc_cells = filter(ExampleAirrTyCHE, location=="germinal_center")
## gcf = formatClones(gc_cells, traits=c("location","sample_time"),
## 	germ="germline_alignment")
## 
## 
## gctrees = getTrees(gcf, build="pml", sub_model="HKY")
## 
## plotTrees(gctrees)[[1]] + geom_tippoint(aes(color=sample_time))


## ----eval=FALSE---------------------------------------------------------------------------------------------------------------
## # edit to your BEAST installation path
## beast <- "/Applications/BEAST 2.7.7/bin/"
## 
## # estimate clock rate of GC B cells
## # if you don't care about convergence, reduce mcmc_length
## # ensure you are providing the correct path to the template file downloaded earlier (see Requirements)
## gctree = getTimeTreesIterate(gcf,
## 	beast=beast,
## 	template="StrictClock_Standard_EmpFreq.xml",
## 	dir="temp",
## 	id="gc_strict",
## 	time="sample_time",
## 	mcmc_length=1e6,
## 	iterations=10,
## 	nproc=2,
## 	CLOCK_RATE_INIT=0.001,
## 	KAPPA_PRIOR_M=0.67,
## 	KAPPA_PRIOR_S=0.2,
## 	ignore=c("freqParameter"))
## 
## 
## gcrate_tree = mean(sapply(gctree$parameters, function(x)filter(x,item=="geneticClockRate")$mean))
## print(gcrate_tree)

## ----echo=FALSE---------------------------------------------------------------------------------------------------------------
print(0.000363)


## ----eval=FALSE---------------------------------------------------------------------------------------------------------------
## gcrate_slope = mean(correlationTest(gctrees, time="sample_time")$slope)
## print(gcrate_slope)


## ----echo=FALSE---------------------------------------------------------------------------------------------------------------
print(0.0003686277)


## ----eval=FALSE---------------------------------------------------------------------------------------------------------------
## mixed_trees <- getTimeTreesIterate(
##   clones,
##   beast    = beast,
##   template = "TraitLinkedExpectedOccupancy_EstTraitClockRates_EmpFreq.xml",
##   trait    = trait,
##   time     = "sample_time",
##   dir      = "temp",
##   id       = "tyche_eo_est",
##   log_every = "auto",
##   nproc     = 2,
##   KAPPA_PRIOR_M = 0.67,
##   KAPPA_PRIOR_S = 0.2,
##   TRAIT_RATE_MEAN_1 = gcrate_tree,
##   TRAIT_RATE_MEAN_2 = 0.000001,
##   TRAIT_RATE_SIGMA_1 = gcrate_tree * 0.01,
##   TRAIT_RATE_SIGMA_2 = 0.001,
##   RATE_INDICATORS = "1 0",
##   TRANSITION_RATE_ALPHA_1 = 0.1,
##   TRANSITION_RATE_ALPHA_2 = 1.0,
##   TRANSITION_RATE_BETA_1  = 0.1,
##   TRANSITION_RATE_BETA_2  = 1.0,
##   log_target   = 2000,
##   mcmc_length  = 1e6,
##   ignore       = c("freqParameter"),
##   iterations   = 20
## )


## ----plot-mixed-trees-result, eval=FALSE--------------------------------------------------------------------------------------
## plotTrees(mixed_trees, scale=10)[[1]] + geom_point(aes(fill=location), pch=21, size=3)


## ----read-beast-output, eval=FALSE--------------------------------------------------------------------------------------------
## mixed_trees <- readBEAST(clones, dir="temp", id="tyche_eo_est", beast=beast, trait=trait)

