## # Enter these commands in a terminal, not an R session!

## 

## # Move to the directory of interest

## mkdir germlines

## 

## # Download the Immcantation repository

## git clone https://bitbucket.org/kleinstein/immcantation

## 

## # Run script to obtain IMGT gapped sequences

## immcantation/scripts/fetch_imgtdb.sh -o germlines

## 

## # View added directories

## ls germlines

## # human  IMGT.yaml  immcantation  mouse  rabbit  rat  rhesus_monkey


## ---- eval=FALSE, warning=FALSE, message=FALSE--------------------------------
## library(dowser)
## library(dplyr)
## 
## data(ExampleAirr)
## 
## # Read in IMGT-gapped sequences
## references = readIMGT(dir = file.path("germlines", "human", "vdj"))
## 
## # remove germline alignment columns for this example
## db = select(ExampleAirr, -"germline_alignment",
##     -"germline_alignment_d_mask")
## 
## # Reconstruct germline sequences
## ExampleAirr = createGermlines(db, references, nproc=1)
## 
## # Check germline of first row
## ExampleAirr$germline_alignment_d_mask[1]
## 
## # "CAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTCAAGCCTGGAGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTC............AGTGACTACTACATGAGCTGGATCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTTCATACATTAGTAGTAGT......AGTAGTTACACAAACTACGCAGACTCTGTGAAG...GGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTGTATTACTGTGCGAGAGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"

