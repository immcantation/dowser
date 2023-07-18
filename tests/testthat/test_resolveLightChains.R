# Load test dataset
#library(testthat)
#library(dowser)
library(dplyr)

data("ExampleMixedDb")

#heavy <- filter(ExampleMixedDb, locus == "IGH")
#light <- filter(ExampleMixedDb, locus != "IGH")

# run resolveLightChains
data <- expect_warning(resolveLightChains(ExampleMixedDb))

expect_equal("seq09" %in% data$sequence_id, FALSE)
expect_equal("seq18" %in% data$sequence_id, FALSE)
expect_equal("seq19" %in% data$sequence_id, FALSE)
expect_equal(data$clone_subgroup[which(data$sequence_id == "seq20")],data$clone_subgroup[which(data$sequence_id == "seq05")])
expect_equal(data$clone_subgroup[which(data$sequence_id == "seq21")],data$clone_subgroup[which(data$sequence_id == "seq15")])
expect_equal(data$clone_subgroup[which(data$sequence_id == "seq26")], 1)
