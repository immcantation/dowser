## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
library(dowser)
library(alakazam)

# load example tsv data
data <- readChangeoDb(file.path("..", "setup", "db_test.tsv"))

# set the clone id 
data$clone_id <- data$expected_clone_id_split_light_T

# seperate the data by heavy and light chains 
heavy <- filter(data, locus == "IGH")
light <- filter(data, locus != "IGH")


# find the clone subgroups 
df <- resolveLightChains(heavy, light, seq = "junction")

print(df$clone_subgroup)
