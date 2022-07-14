# Measurable evolution

**These features are currently only in the development version of Dowser, not the version on CRAN. Please see install page for details on the latest version.**

Dowser implements recently developed phylogenetic tests to detect measurable B cell evolution from longitudinally sampled data.

## Date randomization test

The goal of this test is to determine if a B cell lineage has a detectable relationship between mutation and time. If a lineage is accumulating new mutations over a sample interval, we expect a positive correlation between the divergence (sum of branch length to the most recent common ancestor, MRCA) and time elapsed. 

[Full published deatils on these methods is available here](https://elifesciences.org/articles/70873)

## Set up data structures and trees

This step proceeds as in tree building, but it is important to specify the column of the numeric time point values in the `formatClones` step. In this example we are using simulated data from timepoints 0, 7, and 14 days. Filtering out clones that contain only a single timepoint is not strictly necessary but can improve computing time.

Note it is critical that the clone object is formatted properly using `formatClones` with the `trait` option specified to the sample time. The sample time must be **numeric**.


```r
library(dowser)

# load example AIRR tsv data
data(ExampleAirr)

# Process example data using default settings
clones = formatClones(ExampleAirr, traits="timepoint", minseq=3)

# Column shows timepoints in dataset
print(table(ExampleAirr$timepoint))
#0   7  14 
#62 102 225 

# Calculate number of tissues sampled in tree
timepoints = unlist(lapply(clones$data, function(x)
  length(unique(x@data$timepoint))))

# Filter to multi-type trees
clones = clones[timepoints > 1,]

# Build trees using maximum likelihood (can use alternative builds if desired)
trees = getTrees(clones, build="pml")
```

## Perform date randomization test

Once trees have been built, perform the date randomization test using the function `correlationTest`. By default this test will use a clustered version of the date randomization test that corrects for issues like biased sampling of cell subpopulations, as well as types of sequencing error.



```r
# correlation test with 100 repetitions - in practice use at least 1000
test = correlationTest(trees, permutations=100, time="timepoint")
print(test)

# A tibble: 6 × 12
#  clone_id data    locus  seqs trees    slope     p correlation random_correlat…
#     <dbl> <list>  <chr> <int> <lis>    <dbl> <dbl>       <dbl>            <dbl>
#1     3128 <airrC… N        40 <phy… -0.00205 0.871      -0.173         -0.00953
#2     3184 <airrC… N        12 <phy…  0.00111 0.554       0.649          0.0649 
#3     3140 <airrC… N         9 <phy…  0.00156 0.347       0.630         -0.0278 
#4     3192 <airrC… N         9 <phy…  0.00739 0.564       0.956          0.115  
#5     3115 <airrC… N         6 <phy…  0.00159 0.238       0.565          0.0131 
#6     3139 <airrC… N         6 <phy…  0.00308 0.535       0.821          0.0492 


# use uniform correlaion test (more sensitive, but higher false positive rate)
utest = correlationTest(trees, permutations=100, time="timepoint", perm_type="uniform")
print(utest)

# A tibble: 6 × 12
#  clone_id data   locus  seqs trees    slope      p correlation random_correlat…
#     <dbl> <list> <chr> <int> <lis>    <dbl>  <dbl>       <dbl>            <dbl>
#1     3128 <airr… N        40 <phy… -0.00205 0.832       -0.173         -0.0146 
#2     3184 <airr… N        12 <phy…  0.00111 0.0990       0.649          0.0138 
#3     3140 <airr… N         9 <phy…  0.00156 0.0792       0.630         -0.00741
#4     3192 <airr… N         9 <phy…  0.00739 0.129        0.956          0.00818
#5     3115 <airr… N         6 <phy…  0.00159 0.386        0.565          0.0719 
#6     3139 <airr… N         6 <phy…  0.00308 0.158        0.821         -0.00600
```

The `correlationTest` function returns the same object that was provided as input, but with additional columns. The most important of which is the `p` value column. This gives the p value for the test that the observed correlation between divergence at time (`correlation` column) is greater than in trees with permuted time labels (`random_correlation` shows the mean of randomized correlations). The `slope` column shows the slope of the linear regression line between divergence and time, and represents the mean rate of evolution over time in the lineage. By default, permutations are performed using a clustered permutation procedure detailed in [Hoehn et al. 2021](https://elifesciences.org/articles/70873). This adjusts for confounders like biased sampling at different timepoints and some forms of sequencing error. Uniform permutations, where the timepoint at each tip is permuted independently, can be used by setting `perm_type="uniform"`. This may increase power for small datasets, but will likely increase the false positive rate.
