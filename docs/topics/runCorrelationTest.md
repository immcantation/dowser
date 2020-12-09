**runCorrelationTest** - *Resolve polytomies to have the minimum number of single timepoint clades*

Description
--------------------

`rootToTop` performs root-to-tip regression permutation test


Usage
--------------------
```
runCorrelationTest(
phy,
clone,
permutations,
minlength = 0.001,
polyresolve = TRUE,
permutation = c("clustered", "uniform"),
time = "time",
sequence = "sequence_id",
germline = "Germline",
verbose = TRUE,
alternative = c("greater", "two.sided")
)
```

Arguments
-------------------

phy
:   Tree object

clone
:   airrClone data object corresponding to `phy`

permutations
:   Number of permutations to run

minlength
:   Branch lengths to collapse in trees

polyresolve
:   Resolve polytomies to have a minimum number of 
single timepoint clades

permutation
:   Permute among single timepoint clades or uniformly
among tips

time
:   Column name holding numeric time information

sequence
:   Column name holding sequence ID

germline
:   Germline sequence name

verbose
:   Print lots of rubbish while running?

alternative
:   Is alternative that the randomized correlation are greater than 
or equal to observed, or greater/less than?




Value
-------------------

A list of statistics from running the permutation test.


Details
-------------------

See [correlationTest](correlationTest.md) for details




See also
-------------------

[correlationTest](correlationTest.md).






