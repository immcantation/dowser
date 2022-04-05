# Measurable evolution

Dowser implements recently developed phylogenetic tests to detect measurable B cell evolution from longitudinally sampled data.

## Date randomization test

The goal of this test is to determine if a B cell lineage has a detectable relationship between mutation and time. If a lineage is accumulating new mutations over a sample interval, we expect a positive correlation between the divergence (sum of branch length to the most recent common ancestor, MRCA) and time elapsed. 

[Full published deatils on these methods is available here](https://elifesciences.org/articles/70873)

## Perform date randomization test

This step proceeds as in tree building, but it is important to sprecify the column of the discrete trait you want to analyze in the `formatClones` step. In this example we are using simulated data from nose and lung biospies. However, this could be any discrete trait value such as cell types. Filtering out clones that contain only a single trait value type is not strictly necessary but can dramatically improve computing time.

Note it is critical that the clone object be formatted properly using `formatClones` with the `trait` option specified to the sample time. The sample time must be **numeric**.

By default this test will use a clustered version of the date randomization test that corrects for issues like biased sampling of cell subpopulations, as well as types of sequencing error.


```r
library(dowser)

# Object in which clones have been formatted with numeric timepoint 
# values at the tips, and trees built (see Building Trees section)
data(ExampleClones)

trait = "timepoint"

# get calculate number of timepoint sampled in tree
time_types = unlist(lapply(ExampleClones$data, function(x)
  length(unique(x@data[[trait]]))))

# filter to multi-type trees
ExampleClones = ExampleClones[time_types > 1,]
ExampleClones = ExampleClones[ExampleClones$seqs >= 10,]

# correlation test with 100 repetitions - in practice use at least 1000
test = correlationTest(ExampleClones, permutations=100, time=trait)
print(test)

# A tibble: 2 × 12
#   clone_id data    locus  seqs trees    slope     p correlation random_correlat…
#      <dbl> <list>  <chr> <int> <lis>    <dbl> <dbl>       <dbl>            <dbl>
# 1     3128 <airrC… N        43 <phy… -0.00146 0.693      -0.122          -0.0237
# 2     3184 <airrC… N        12 <phy…  0.00115 0.475       0.748          -0.0449


# use uniform correlaion test (more sensitive, but higher false positive rate)
utest = correlationTest(ExampleClones, permutations=100, time=trait, permutation="uniform")
print(utest)

# A tibble: 2 × 12
#   clone_id data    locus  seqs trees    slope     p correlation random_correlat…
#      <dbl> <list>  <chr> <int> <lis>    <dbl> <dbl>       <dbl>            <dbl>
# 1     3128 <airrC… N        43 <phy… -0.00146 0.861      -0.122          0.00934
# 2     3184 <airrC… N        12 <phy…  0.00115 0.109       0.748         -0.0440 
```


As shown above, the `correlationTest` function returns the same object that was provided as input, but with additional columns. The most important of which is the `p` value column. This gives the p value for the test that the observed correlation between divergence at time (`correlation` column) is greater than in trees with permuted time labels (`random_correlation` shows the mean of randomized correlations. The `slope` column shows the slope of the linear regression line between divergence and time, and represents the mean rate of evolution over time in the lineage. By default, permutations are performed using a clustered permutation procedure detailed in [Hoehn et al. 2021](https://elifesciences.org/articles/70873). This adjusts for confounders like biased sampling at different timepoints and some forms of sequencing error. Uniform permutations, where the timepoint at each tip is permuted independently, can be used by setting `permutaion="uniform"`. This may increase power for small datasets, but will likely increase the false positive rate.
