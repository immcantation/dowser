## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
library(dowser)

# load example AIRR tsv data
data(ExampleAirr)

# Subset data for this example
ExampleAirr = ExampleAirr[ExampleAirr$clone_id %in% c("3170", "3184"),]
ExampleAirr$subject_id = "Subject_1"

# Process example data using default settings
clones = formatClones(ExampleAirr)

print(clones)

# Process example data keeping samples from different times
# distinct, adding duplicate_count among collapsed sequences,
# and show the sample_id within each clone in the tibble.
clones = formatClones(ExampleAirr, traits=c("sample_id","c_call"),
    num_fields=c("duplicate_count"), columns=c("subject_id"))

print(clones)

## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
#  
#  clones = getTrees(clones, nproc=1)
#  
#  print(clones)
#  ## A tibble: 2 x 6
#  #  clone_id data       locus  seqs subject_id trees
#  #     <dbl> <list>     <chr> <int> <chr>      <list>
#  #1     3170 <airrClon> N        13 Subject_1  <phylo>
#  #2     3184 <airrClon> N        12 Subject_1  <phylo>

## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
#  # exec here is set to dnapars position in the Docker image.
#  clones = getTrees(clones, build="dnapars", exec="/usr/local/bin/dnapars", nproc=1)
#  
#  clones
#  ## A tibble: 2 x 6
#  #  clone_id data       locus  seqs subject_id trees
#  #     <dbl> <list>     <chr> <int> <chr>      <list>
#  #1     3170 <airrClon> N        13 Subject_1  <phylo>
#  #2     3184 <airrClon> N        12 Subject_1  <phylo>

## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
#  
#  clones = getTrees(clones, build="pml")
#  
#  print(clones)
#  ## A tibble: 2 x 6
#  #  clone_id data       locus  seqs subject_id trees
#  #     <dbl> <list>     <chr> <int> <chr>      <list>
#  #1     3170 <airrClon> N        13 Subject_1  <phylo>
#  #2     3184 <airrClon> N        12 Subject_1  <phylo>

## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
#  
#  # exec here is set to dnaml position in the Docker image.
#  clones = getTrees(clones, build="dnaml", exec="/usr/local/bin/dnaml")
#  
#  clones
#  # A tibble: 2 x 6
#  #  clone_id data       locus  seqs subject_id trees
#  #     <dbl> <list>     <chr> <int> <chr>      <list>
#  #1     3170 <airrClon> N        13 Subject_1  <phylo>
#  #2     3184 <airrClon> N        12 Subject_1  <phylo>

## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
#  
#  # exec here is set to raxml position in the Docker image.
#  clones = getTrees(clones, build="raxml", exec="/usr/local/bin/raxml-ng")
#  
#  clones
#  # A tibble: 2 x 6
#  #  clone_id data       locus  seqs subject_id trees
#  #     <dbl> <list>     <chr> <int> <chr>      <list>
#  #1     3170 <airrClon> N        13 Subject_1  <phylo>
#  #2     3184 <airrClon> N        12 Subject_1  <phylo>

## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
#  
#  # exec here is set to IgPhyML position in the Docker image.
#  clones = getTrees(clones, build="igphyml",
#      exec="/usr/local/share/igphyml/src/igphyml", nproc=1)
#  
#  print(clones)
#  ## A tibble: 2 x 7
#  #  clone_id data       locus  seqs subject_id trees        parameters
#  #     <dbl> <list>     <chr> <int> <chr>      <named list> <named list>
#  #1     3170 <airrClon> N        13 Subject_1  <phylo>      <named list [13]>
#  #2     3184 <airrClon> N        12 Subject_1  <phylo>      <named list [13]>
#  
#  
#  clones$parameters[[1]]$omega_mle
#  #[1] 0.5286

## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
#  
#  # exec here is set to IgPhyML position in the Docker image.
#  clones = getTrees(clones, build="igphyml",
#      exec="/usr/local/share/igphyml/src/igphyml", nproc=1, partition="cf")
#  
#  print(clones)
#  ## A tibble: 2 x 7
#  #  clone_id data       locus  seqs subject_id trees        parameters
#  #     <dbl> <list>     <chr> <int> <chr>      <named list> <named list>
#  #1     3170 <airrClon> N        13 Subject_1  <phylo>      <named list [13]>
#  #2     3184 <airrClon> N        12 Subject_1  <phylo>      <named list [13]>

