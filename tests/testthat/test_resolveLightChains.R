# Load test dataset
library(alakazam)
#load(file.path("test", "setup", "db_test.tsv"), envir=e1)
data <- readChangeoDb(file.path("..", "setup", "db_test.tsv"))

data$clone_id <- data$expected_clone_id_split_light_T
heavy <- filter(data, locus == "IGH")
light <- filter(data, locus != "IGH")

# run resolveLightChains
data <- resolveLightChains(heavy, light, seq = "junction")

expect_equal("seq09" %in% data$sequence_id, FALSE)
expect_equal("seq18" %in% data$sequence_id, FALSE)
expect_equal("seq19" %in% data$sequence_id, FALSE)
expect_equal(data$clone_subgroup[which(data$sequence_id == "seq20")],data$clone_subgroup[which(data$sequence_id == "seq05")])
expect_equal(data$clone_subgroup[which(data$sequence_id == "seq21")],data$clone_subgroup[which(data$sequence_id == "seq15")])
expect_equal(data$clone_subgroup[which(data$sequence_id == "seq26")], 1)
