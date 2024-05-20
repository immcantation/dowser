#library(testthat)
#library(dowser)

# Test the correlationTest
test_that("correlationTest", {
  data(ExampleAirr)
  db = formatClones(ExampleAirr, traits="timepoint", minseq=3)
  timepoints = unlist(lapply(db$data, function(x)
    length(unique(x@data$timepoint))))
  
  # Filter to multi-type trees
  db = db[timepoints > 1,]
  
  # Build trees using maximum likelihood (can use alternative builds if desired)
  trees = getTrees(db[3:6,], build="pml")
  
  # do the correlationTest
  trees <- correlationTest(trees, time = "timepoint")
  expect_equal(trees$nclust,
               c(3,2,4,2))
  expect_equal(trees$clone_id,
               c(3140, 3192, 3115, 3139))
})


