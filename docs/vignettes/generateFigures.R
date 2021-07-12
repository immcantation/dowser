library(dowser)
library(ggtree)
library(ggplot2)
library(dplyr)
library(gridExtra)

data(ExampleDb)

clones = formatClones(ExampleDb, traits=c("c_call"),
    num_fields=c("duplicate_count"), columns=c("sample_id"),
    minseq=10)

# build maximum parsimony trees
clones = getTrees(clones)

plots = plotTrees(clones)

#Plot the largest tree
plots[[1]]

png("figures/Plotting-Vignette-basic.png")
plots[[1]]
dev.off()

# Re-scale branches to represent mutations between nodes
clones_mutations = scaleBranches(clones, edge_type="mutations")

# Plot, set scale bar to represent 10 mutations
plots = plotTrees(clones_mutations, scale=10)

#Plot the largest tree
plots[[1]]

png("figures/Plotting-Vignette-scaled.png")
plots[[1]]
dev.off() 

# Plot tree with sequence isotype at the tips.
plots = plotTrees(clones, tips="c_call")

#Plot the largest tree
plots[[1]]

png("figures/Plotting-Vignette-c_call.png")
plots[[1]]
dev.off()

# Plot tree with sequence isotype at the tips, with sizes set to number of duplicates
plots = plotTrees(clones, tips="c_call", tipsize="duplicate_count")

#Plot the largest tree
plots[[1]]

png("figures/Plotting-Vignette-c_call_duplicate.png")
plots[[1]]
dev.off()

# These calls create the same plot:

# Plot tree with sequence isotype at the tips, with palette "Set1"
plots = plotTrees(clones, tips="c_call", tipsize=2,
    tip_palette="Set1")

# or, specify a named palette vector
custom_palette=c("IGHA"="#E41A1C", "IGHG"="#377EB8",
    "IGHD"="#4DAF4A", "Germline"="#984EA3")
plots = plotTrees(clones, tips="c_call", tipsize=2,
    tip_palette=custom_palette)

# or, use the getPalette function to create a named palette vector
custom_palette=getPalette(c("IGHA", "IGHG", "IGHD", "Germline"), "Set1")
plots = plotTrees(clones, tips="c_call", tipsize=2,
    tip_palette=custom_palette)

#Plot the largest tree
plots[[1]]

png("figures/Plotting-Vignette-c_call_set1.png")
plots[[1]]
dev.off()

plots = plotTrees(clones, tips="c_call", tipsize=2)

#Plot the largest tree
treeplot = plots[[1]] + geom_tiplab() + 
    geom_vline(xintercept=c(0.05,0.1,0.15,0.2,0.25),
        linetype="dashed",color="grey") + xlim(0,0.3) +
    ggtitle("Example B cell tree")

treeplot

png("figures/Plotting-Vignette-ggtree.png")
treeplot
dev.off()

plots = plotTrees(clones, tips="c_call", tipsize=2)

# pass any arguments you would pass to grDevices::pdf
treesToPDF(plots, file="trees.pdf", nrow=2, ncol=2)

# pass any arguments you would pass to grDevices::pdf
png("figures/Plotting-Vignette-all.png")
grid.arrange(grobs=plots[1:4],ncol=2)
dev.off()


clones = collapseNodes(clones)
plots = plotTrees(clones, tips="c_call", tipsize=2,
    node_nums=TRUE, labelsize=7)

#treesToPDF(plots, file="trees.pdf", nrow=2, ncol=2)

# pass any arguments you would pass to grDevices::pdf
png("figures/Sequences-Vignette-all.png")
grid.arrange(grobs=plots[1],ncol=1)
dev.off()

getSeq(clones, node=54, clone=3128)


# cd <data directory>
# 
# # Download the Immcantation repository
# git clone bitbucket.org/kleinstein/immcantation
# 
# # Run script to obtain IMGT gapped sequences
# Immcantation/scripts/fetch_imgtdb.sh

library(dowser)
data(ExampleDb)

references = readIMGT(dir = "germlines/human/vdj")

# remove germline alignment columns for this example
ExampleDb = select(ExampleDb, -"germline_alignment", 
    -"germline_alignment_d_mask")

# Reconstruct germline sequences
ExampleDb = createGermlines(ExampleDb,references,nproc=1)

# Check germline of first row
ExampleDb$germline_alignment_d_mask[1]
