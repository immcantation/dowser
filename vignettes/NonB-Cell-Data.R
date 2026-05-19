## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
library(dowser)
data(ExampleAirr)
ExampleAirr <- dplyr::select(dplyr::filter(ExampleAirr, clone_id=="3128"),
                             sequence_id, sequence_alignment, germline_alignment)
clones <- formatClones(ExampleAirr, germ="germline_alignment")
trees <- getTrees(clones)

