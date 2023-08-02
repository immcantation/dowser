# Load test dataset
#library(testthat)
#library(dowser)

data("ExampleMixedDb")

# run resolveLightChains
data <- resolveLightChains(ExampleMixedDb)

expect_equal(data$clone_subgroup[which(data$sequence_id == "d12_BCell_GTATTCTTCACTGGGC-1")], 2)
expect_equal(data$clone_subgroup[which(data$sequence_id == "d12_FNA_CGATTGACAAGGTTCT-1")],
             data$clone_subgroup[which(data$sequence_id == "d60_FNA_TACACGACATGACATC-1")])
expect_equal(data$clone_subgroup[which(data$sequence_id == "d12_BCell_GTATTCTTCACTGGGC-1")],
             data$clone_subgroup[which(data$sequence_id == "d12_BCell@GTATTCTTCACTGGGC-1_contig_1")])
expect_equal(data$clone_subgroup[which(data$sequence_id == "d60_FNA_CTCGGGATCTGGCGTG-1")],
             data$clone_subgroup[which(data$sequence_id == "d60_FNA@CTCGGGATCTGGCGTG-1_contig_1")])
expect_equal(data$clone_subgroup_id[which(data$sequence_id == "d60_FNA_CTCGGGATCTGGCGTG-1")],
             "45155_2")
expect_equal(data$clone_subgroup[which(data$sequence_id == "d60_BCell_ACTATCTAGACCTAGG-1")],
             as.integer(strsplit(data$clone_subgroup_id[which(data$sequence_id == "d60_BCell_ACTATCTAGACCTAGG-1")], "_")[[1]][2]))



