**rootToTip** - *Performs root-to-tip regression test on set of trees*

Description
--------------------

`rootToTop` performs root-to-tip regression permutation test


Usage
--------------------
```
rootToTip(
trees,
time = "time",
permutations = 1000,
germline = "Germline",
minlength = 0.001,
alternative = c("two.sided", "greater", "less")
)
```

Arguments
-------------------

trees
:   Tibble with trees and data

time
:   Column id for timepoint

permutations
:   Number of permutations for test

germline
:   Germline sequence name

minlength
:   Branch lengths to collapse in trees

alternative
:   Perform one-sided (`greater` or `less`)
or `two.sided` test




Value
-------------------

A `tibble` with pearson correlation between divergene
and time, mean permuted correlation, p value(s), number of permutations,
and number of sequences


Details
-------------------

Output data table columns:
clone_id = clone id
observed = observed pearson correlation
permuted = mean permuted correlation
pgt = p value for DELTA < 0
plt = p value for DELTA > 0




See also
-------------------

Uses output from [getTrees](getTrees.md).






