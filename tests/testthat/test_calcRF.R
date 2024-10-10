# Load test dataset
#library(testthat)
#library(dowser)
library(dplyr)

data("ExampleClones")

# read in the phangorn:pratchet tree
tree1 <- ExampleClones$trees[[which(ExampleClones$clone_id == 3157)]]

# get the phangorn::pml tree -- should be the same 
tree2 <- getTrees(filter(ExampleClones, clone_id == 3157), build = "pml")
tree2 <- tree2$trees[[1]]

# test calcRF
rf_dist <- calcRF(tree1, tree2)

# test the seq names in both trees
# test the dist 
expect_equal(unique(sort(tree1$tip.label) == sort(tree2$tip.label)), TRUE)
expect_equal(rf_dist, 0)

# now do one where you would know the difference -- should be 2
data <- ExampleClones[28,]
data <- getTrees(data)
data2 <- data
data2$data[[1]]@data$sequence[which(data2$data[[1]]@data$sequence_id == "GN5SHBT06HS8XB")] <- 
  data2$data[[1]]@data$sequence[which(data2$data[[1]]@data$sequence_id == "GN5SHBT01DGOTM")]
data2 <- getTrees(data2)

# get the distance 
rf_dist <- calcRF(data$trees[[1]], data2$trees[[1]])

# test the tip names and the distance expectation 
expect_equal(unique(sort(data$trees[[1]]$tip.label) == sort(data2$trees[[1]]$tip.label)), TRUE)
expect_equal(rf_dist, 2)