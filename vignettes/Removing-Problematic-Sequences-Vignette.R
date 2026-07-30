library(dowser)
library(dplyr)
library(ggtree)

data(ExampleAirr)

# Simulate a partial/truncated sequence
row <- which(ExampleAirr$sequence_id == "GN5SHBT04BIIPE")
orig <- ExampleAirr$sequence_alignment[row]
half <- floor(nchar(orig)/2)
ExampleAirr$sequence_alignment[row] <- paste0(substr(orig, 1, half),
  paste(rep("N", nchar(orig) - half), collapse=""))

filtered <- filterPartialSeqs(ExampleAirr)

"GN5SHBT04BIIPE" %in% filtered$sequence_id

data(ExampleAirr)

# Simulate an oversequencing/PCR error artifact within clone 3128
ExampleAirr$duplicate_count[ExampleAirr$sequence_id == "GN5SHBT07H5AOD"] = 1000
ExampleAirr$duplicate_count[ExampleAirr$sequence_id == "GN5SHBT08H9MGK"] = 101

clone <- filter(ExampleAirr, clone_id == 3128)
clones <- formatClones(clone, trait="biopsy", num_fields="duplicate_count")
trees <- scaleBranches(getTrees(clones))

plotTrees(trees, tipsize="duplicate_count", tips="biopsy", scale=1)[[1]] +
  geom_tiplab(size=3) + xlim(0, 55)

filtered <- filterCombs(ExampleAirr, dup_count_thresh=100, exponent=1)

fclone <- filter(filtered, clone_id == 3128)
fclones <- formatClones(fclone, trait="biopsy", num_fields="duplicate_count")
ftrees <- scaleBranches(getTrees(fclones))

plotTrees(ftrees, tipsize="duplicate_count", tips="biopsy", scale=1)[[1]] +
  geom_tiplab(size=3) + xlim(0, 55)
