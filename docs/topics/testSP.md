**testSP** - *Performs SP (switch proportion) test on switch data*

Description
--------------------

`testSP` performs an SP test


Usage
--------------------
```
testSP(
switches,
permuteAll = FALSE,
from = NULL,
to = NULL,
dropzeros = TRUE,
bylineage = FALSE,
pseudocount = 0,
alternative = c("two.sided", "greater", "less"),
tip_switch = 20,
exclude = FALSE
)
```

Arguments
-------------------

switches
:   Data frame from bootstrapTrees

permuteAll
:   Permute among trees?

from
:   Include only switches from this state?

to
:   Include only switches to this state?

dropzeros
:   Drop switches with zero counts?

bylineage
:   Perform test for each lineage individually?

pseudocount
:   Pseudocount for P value calculations

alternative
:   Perform one-sided (`greater` or `less`)
or `two.sided` test

tip_switch
:   maximum tip/switch ratio

exclude
:   exclude clones with tip/switch ratio > `tip_switch`?




Value
-------------------

A list containing a `tibble` with mean SP statistics, and another 
with SP statistics per repetition.


Details
-------------------

Output data table columns:
RECON = SP for observed data
PERMUTE = SP for permuted data
DELTA = RECON - PERMUTE
PLT = p value for DELTA < 0
PGT = p value for DELTA < 0

+ `FROM`: State going from.
+ `TO`: State going to.
+ `RECON`: SP for observed data.
+ `PERMUTE`: SP for permuted data.
+ `DELTA`:  RECON - PERMUTE.
+ `PLT`: p value that DELTA < 0
+ `PGT`: p value that DELTA > 0
+ `STAT`: Statistic used (SP).
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
# btrees = bootstrapTrees(clones[1:2], bootstraps=10, nproc=1,
# igphyml=igphyml, trait="sample_id")
# testSP(btrees$switches)
```



See also
-------------------

Uses output from [bootstrapTrees](bootstrapTrees.md). Related to [testPS](testPS.md)
and [testSC](testSC.md).






