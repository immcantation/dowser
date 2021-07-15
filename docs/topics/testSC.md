**testSC** - *Performs SC (switch count) test on switch data*

Description
--------------------

`testSC` performs an SC test


Usage
--------------------
```
testSC(
switches,
dropzeros = TRUE,
bylineage = FALSE,
pseudocount = 0,
from = NULL,
to = NULL,
permuteAll = FALSE,
alternative = c("two.sided", "greater", "less")
)
```

Arguments
-------------------

switches
:   Data frame from bootstrapTrees

dropzeros
:   Drop switches with zero counts?

bylineage
:   Perform test for each lineage individually?

pseudocount
:   Pseudocount for P value calculations

from
:   Include only switches from this state?

to
:   Include only switches to this state?

permuteAll
:   Permute among trees?

alternative
:   Perform one-sided (`greater` or `less`)
or `two.sided` test




Value
-------------------

A list containing a `tibble` with mean SC statistics, and another 
with SC statistics per repetition.


Details
-------------------

Output data table columns:
RECON = SC for observed data
PERMUTE = SC for permuted data
DELTA = RECON - PERMUTE
PLT = p value for DELTA < 0
PGT = p value for DELTA < 0

+ `FROM`: State going from.
+ `TO`: State going to.
+ `RECON`: SC for observed data.
+ `PERMUTE`: SC for permuted data.
+ `DELTA`:  RECON - PERMUTE.
+ `PLT`: p value that DELTA < 0
+ `PGT`: p value that DELTA > 0
+ `STAT`: Statistic used (SC).
+ `REP`: Bootstrap repetition.
+ `REPS`: Total number of ootstrap repetition.




Examples
-------------------

```R
### Not run:
igphyml <- "~/apps/igphyml/src/igphyml"
# data(ExampleAirr)
# ExampleAirr$sample_id = sample(ExampleAirr$sample_id)
# clones = formatClones(ExampleAirr, trait="sample_id")
# btrees = bootstrapTrees(clones[1:2], bootstraps=100, nproc=1,
# igphyml=igphyml, trait="sample_id", id="temp", dir="temp")
# testSC(btrees$switches)
```



See also
-------------------

Uses output from [bootstrapTrees](bootstrapTrees.md). Related to [testPS](testPS.md)
and [testSP](testSP.md).






