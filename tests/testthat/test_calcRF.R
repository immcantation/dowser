# Load test dataset
#library(testthat)
#library(dowser)
library(dplyr)

data("ExampleClones")

# read in the phangorn:pratchet tree
tree1 <- ExampleClones$trees[[which(ExampleClones$clone_id == 3157)]]

# get the phangorn::pml tree
tree2 <- getTrees(filter(ExampleClones, clone_id == 3157), build = "pml")
tree2 <- tree2$trees[[1]]

# test calcRF
rf_dist <- calcRF(tree1, tree2)

# test the seq names in both trees
# test the dist 
expect_equal(unique(sort(tree1$tip.label) == sort(tree2$tip.label)), TRUE)
expect_equal(rf_dist, 2)