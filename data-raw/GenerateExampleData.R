# Generate example trees

# Imports
library(dowser)
library(alakazam)
library(dplyr)
library(igraph)
library(Matrix)

set.seed(1)

#### Helper functions

#draw the next state based on a weighted substitution matrix
simStateQ = function(pstate,length,rate,Q){
    P = expm(Q*length*rate)
    state = sample(1:nrow(Q),size=1,prob=P[pstate,])
    return(state)
}

#simulate metadata along a tree based on a Q matrix
#Q = matrix(0,nrow=2,ncol=2)
#Q[1,2] = 1
#Q[2,1] = 2
#diag(Q) = -rowSums(Q)
#m = sum(freqs * -diag(Q))
#Q = Q/m
simulateSubs = function(tree,rate,freqs,Q,states,start=NULL){
    if(is.null(start)){
        start = sample(1:length(states),size=1,prob=freqs)
    }
    edges = tree$edge
    treestates = rep(-1,length=max(edges))
    uca = ape::getMRCA(tree,tip=tree$tip.label )
    parent_nodes = uca
    treestates[uca] = start
    nchild = length(edges[edges[,1] == parent_nodes[1],2])
    while(nchild > 0){ #stolen from shazam::shmulateTree
        new_parents = c()
        for(p in parent_nodes){
            children = edges[edges[,1] == p,2]
            pstate = treestates[p]
            for(ch in children){
                new_parents = union(new_parents, ch)
                edge = which(edges[,1] == p & edges[,2] == ch)
                length = tree$edge.length[edge]
                treestates[ch] = simStateQ(pstate,length,rate,Q)
            }
        }
        parent_nodes = new_parents
        nchild = length(edges[edges[,1] %in% parent_nodes,2])
    }
    tree$state = states[treestates]
    tree$node.label = tree$state[(length(tree$tip.label)+1):length(tree$state)]
    return(tree)
}


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

# Simulate discrete traits
tistates = c("Nose", "Lung")
freqs = c(1,0)
Q = matrix(0,nrow=2,ncol=2)
Q[1,2] = 10
Q[2,1] = 1
diag(Q) = -rowSums(Q)
m = sum(freqs * -diag(Q))
Q = Q/m

simtrees = lapply(trees$trees,function(x) 
        x=simulateSubs(x,10,freqs,Q,tistates))

simulated = data.frame()
for(tree in simtrees){
    s2l = tree$state[1:length(tree$tip.label)]
    names(s2l) = tree$tip.label
    tclone = ExampleAirr[ExampleAirr$sequence_id %in% 
        tree$tip.label & ExampleAirr$clone_id == tree$name,]
    tclone$location = s2l[tclone$sequence_id]
    simulated = rbind(simulated,tclone)
}
simulated$biopsy = simulated$location

# Simulate timepoints
tistates = c("0d", "7d", "14d")
freqs = c(1,0,0)
Q = matrix(0,nrow=3,ncol=3)
Q[1,2] = 1
Q[2,3] = 1
diag(Q) = -rowSums(Q)
m = sum(freqs * -diag(Q))
Q = Q/m

simtrees = lapply(trees$trees,function(x) 
        x=simulateSubs(x,30,freqs,Q,tistates))

simulated2 = data.frame()
for(tree in simtrees){
    s2l = tree$state[1:length(tree$tip.label)]
    names(s2l) = tree$tip.label
    tclone = simulated[simulated$sequence_id %in% 
        tree$tip.label & simulated$clone_id == tree$name,]
    tclone$location = s2l[tclone$sequence_id]
    simulated2 = rbind(simulated2,tclone)
}
simulated2$timepoint = as.numeric(gsub("d","",simulated2$location))


ExampleAirr = simulated2

f = formatClones(ExampleAirr, 
    traits=c("c_call","biopsy","timepoint"), 
    num_fields="duplicate_count")

ExampleClones <- getTrees(f)

# Save
usethis::use_data(ExampleAirr, overwrite=TRUE)
usethis::use_data(ExampleClones, overwrite=TRUE)

