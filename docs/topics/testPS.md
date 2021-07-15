**testPS** - *Performs PS (parsimony score) test on switch data*

Description
--------------------

`testPS` performs a PS test


Usage
--------------------
```
testPS(
switches,
bylineage = FALSE,
pseudocount = 0,
alternative = c("less", "two.sided", "greater")
)
```

Arguments
-------------------

switches
:   Data frame from bootstrapTrees

bylineage
:   Perform test for each lineage individually? (FALSE)

pseudocount
:   Pseudocount for P value calculations

alternative
:   Perform one-sided (`greater` or `less`)
or `two.sided` test




Value
-------------------

A list containing a `tibble` with mean PS statistics, and another 
with PS statistics per repetition.


Details
-------------------

Output data table columns:
RECON = PS for observed data
PERMUTE = PS for permuted data
DELTA = RECON - PERMUTE
PLT = p value for DELTA < 0
PGT = p value for DELTA < 0

+ `RECON`: PS for observed data.
+ `PERMUTE`: PS for permuted data.
+ `DELTA`:  RECON - PERMUTE.
+ `PLT`: p value that DELTA < 0
+ `PGT`: p value that DELTA > 0
+ `STAT`: Statistic used (PS).
+ `REP`: Bootstrap repetition.
+ `REPS`: Total number of ootstrap repetition.




Examples
-------------------

```R
### Not run:
igphyml <- "~/apps/igphyml/src/igphyml"
# data(ExampleAirr)
# ExampleAirr$sample_id <- sample(ExampleAirr$sample_id)
# clones <- formatClones(ExampleAirr, trait="sample_id")
# btrees <- bootstrapTrees(clones[1:2], bootstraps=10, nproc=1,
# igphyml=igphyml, trait="sample_id")
# testPS(btrees$switches)
```



See also
-------------------

Uses output from [bootstrapTrees](bootstrapTrees.md). Related to [testSP](testSP.md)
and [testSC](testSC.md).






