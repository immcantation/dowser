## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
library(dowser)

data("ExampleMixedDb")

# find the clone subgroups 
ExampleMixedDb <- resolveLightChains(ExampleMixedDb)
print(ExampleMixedDb$clone_subgroup)

# since the example data has a "bulk" heavy chain sequence, the missing cell_id will need to be replaced
ExampleMixedDb$cell_id[16] <- "bulk"
# for paired data, ensure that the chain is HL. collapse = F here to get multiple clones due to the sequences being very similiar. 
clones <- formatClones(ExampleMixedDb, chain="HL", nproc=1, collapse = F)
print(clones)

# this is build with the default tree building method in dowser
# phangorn's pratchet method
clones <- getTrees(clones, nproc = 1)

print(clones)