# Generate example trees

# Imports
library(dowser)
library(alakazam)
library(dplyr)
library(igraph)

#### Generate example database ####

# Load data
ExampleAirr <- readChangeoDb("data-raw/ExampleAIRR.tsv")
ExampleAirr2 <- readChangeoDb("data-raw/ExampleDb.gz")

m = match(ExampleAirr$sequence_id, ExampleAirr2$sequence_id)
mean(ExampleAirr$sequence_id == ExampleAirr$sequence_id[m])
ExampleAirr$germline_alignment_d_mask = ExampleAirr2$germline_alignment_d_mask[m]
ExampleAirr$v_call_genotyped = ExampleAirr2$v_call_genotyped[m]
ExampleAirr$duplicate_count = ExampleAirr2$duplicate_count[m]
ExampleAirr$clone_id = ExampleAirr2$clone_id[m]
ExampleAirr$sample_id = ExampleAirr2$sample_id[m]


ExampleAirr <- ExampleAirr[c("sequence_id",
                         "sequence_alignment",
                         "germline_alignment",
                         "germline_alignment_d_mask",
                         "rev_comp",
                         "productive",
                         "v_call",
                         "v_call_genotyped",
                         "d_call",
                         "j_call",
                         "c_call",
                         "junction",
                         "junction_length",
                         "np1_length",
                         "np2_length",
                         "duplicate_count",
                         "clone_id",
                         "sample_id",
                         "v_germline_start",
                         "v_germline_end",
                         "d_germline_start",
                         "d_germline_end",
                         "j_germline_start",
                         "j_germline_end"
                         )]

c_trans <- c(IGHM="IgM", IGHD="IgD", IGHA="IgA", IGHG="IgG")
ExampleAirr <- ExampleAirr %>%
    mutate(c_call=translateStrings(c_call, c_trans),
           germline_alignment=germline_alignment_d_mask)

refs = readIMGT("~/share/germlines/imgt/human/vdj")

ExampleAirr = createGermlines(ExampleAirr, refs)

clones = formatClones(ExampleAirr, trait="c_call", num_fields="duplicate_count")

# Build maxmimum parsimony trees for first two clones using 
# phangorn package in R
trees <- getTrees(clones)

ExampleClones = trees

# Save
usethis::use_data(ExampleAirr, overwrite=TRUE)
usethis::use_data(ExampleClones, overwrite=TRUE)
