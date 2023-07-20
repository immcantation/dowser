**correlationTest** - *Run date randomization test for temporal signal on a set of trees.*

Description
--------------------

`correlationTest` performs root-to-tip regression date randomization test


Usage
--------------------
```
correlationTest(
clones,
permutations = 1000,
minlength = 0.001,
perm_type = c("clustered", "uniform"),
time = "time",
sequence = "sequence_id",
germline = "Germline",
verbose = FALSE,
polyresolve = TRUE,
alternative = c("greater", "two.sided"),
storeTree = FALSE,
nproc = 1
)
```

Arguments
-------------------

clones
:   A `tibble` object containing airrClone and `phylo` objects

permutations
:   Number of permutations to run

minlength
:   Branch lengths to collapse in trees

perm_type
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

polyresolve
:   Resolve polytomies to have a minimum number of 
single timepoint clades

alternative
:   Is alternative that the randomized correlation are greater than 
or equal to observed, or greater/less than?

storeTree
:   Store the tree used?

nproc
:   Number of cores to use for calculations. Parallelizes by tree.




Value
-------------------

A `tibble` with the same columns as clones, but additional
columns corresponding to test statistics for each clone.


Details
-------------------

Object returned contains these columns which are added or modified from input:
 
+ `data`: airrClone object, same as input but with additional columns 
"cluster" which correspond to permutation cluster, and "divergence."
+ `slope`: Slope of linear regression between divergence and time.
+ `correlation`: Correlation between divergence and time.
+ `p`: p value of correlation compared to permuted correlations.
+ `random_correlation`: Mean correlation of permutation replicates.
+ `min_p`: Minimum p value of data, determined by either the number of
distinct clade/timepoint combinations or number of permutations.
+ `nposs`: Number of possible distinct timepoint/clade combinations.
+ `nclust`: Number of clusters used in permutation. If perm_type="uniform"
this is the number of tips.
+ `p_gt/p_lt`: P value that permuted correlations are greater or less 
than observed correlation. Only returned if alternative = "two.sided"
+ `test_trees`:  The [phylo](http://www.rdocumentation.org/packages/ape/topics/read.tree) tree objects used, possibly with
resolved polytomies.





See also
-------------------

Uses output from `getTrees`.






