# Functions for performing discrete trait analysis on B cell lineage trees

# Write a clone's sequence alignment to a fasta file
# 
# \code{cloneToFasta} write clone sequences as a fasta file
# @param    c            airrClone object
# @param    fastafile    file to be exported
# @param    germid       sequence id of germline
# @param    trait        trait to include in sequence ids
# @param    empty        don't include real sequence information
#
# @return   Name of exported fasta file.
cloneToFasta <- function(c, fastafile, germid, trait=NULL, empty=FALSE){
  text <- ""
  if(!is.null(trait)){
    c@data$sequence_id <- paste(c@data$sequence_id,c@data[,trait],sep="_")
  }
  for(i in 1:nrow(c@data)){
    text <- paste0(text,">",c@data[i,]$sequence_id,"\n")
    if(!empty){
      if(c@phylo_seq == "sequence"){
        text <- paste0(text,c@data[i,]$sequence,"\n")
      }else if(c@phylo_seq == "lsequence"){
        text <- paste0(text,c@data[i,]$lsequence,"\n")
      }else if(c@phylo_seq == "hlsequence"){
        text <- paste0(text,c@data[i,]$hlsequence,"\n")
      }else{
        stop(paste("phylo_seq not recognized",c@clone))
      }
    }else{
      text <- paste0(text,"ATG\n")
    }
  }
  text <- paste0(text,">",germid,"\n")
  if(!empty){
    if(c@phylo_seq == "sequence"){
      text <- paste0(text,c@germline,"\n")
    }else if(c@phylo_seq == "lsequence"){
      text <- paste0(text,c@lgermline,"\n")
    }else if(c@phylo_seq == "hlsequence"){
      text <- paste0(text,c@hlgermline,"\n")
    }else{
      stop(paste("phylo_seq not recognized",c@clone))
    }
  }else{
    text <- paste0(text,"ATG\n")
  }
  write(text,file=fastafile,append=FALSE)
  return(fastafile)
}


#' Read a fasta file into a list of sequences
#' \code{readFasta} reads a fasta file
#' @param    file      FASTA file
#'
#' @return   List of sequences
#' @export
readFasta <- function(file){
  f <- readLines(file)
  if(length(f) == 1){
    return(NULL)
  }
  seqs <- list()
  id <- NA
  for(line in f){
    if(grepl("^>",line)){
      id <- gsub(">","",line)
      seqs[[id]] <- ""
    }else{
      if(is.na(id)){
        stop(paste("Error reading",file))
      }
      seqs[[id]] <- paste0(seqs[[id]],line)
    }
  }
  seqs
}

#' Read in a parsimony model file
#' 
#' \code{readModelFile} Filler
#' @param    file         parimony model file.
#' @param    useambig     use ambiguous naming as specified in the file?
#'
#' @return   A named vector containing the states of the model
#'
#' @seealso \link{makeModelFile}, \link{findSwitches}, \link{getTrees}
#'
#' @export
readModelFile <- function(file, useambig=FALSE){
  #set up pallete
  mfile <- readLines(file)
  mstart <- which(mfile == "#STATES")
  mend <- which(mfile == "")
  mend <- min(mend[which(mfile == "") > mstart])
  states <- mfile[(mstart+1):(mend-1)]
  names(states) <- states
  
  if(useambig){
    stop("Ambiguous states not currently supported, please useambig=FALSE")
    astart <- which(mfile == "#AMBIGUOUS")
    ambigs <- mfile[(astart+1):length(mfile)]
    asplit <- strsplit(ambigs, split=" ")
    ambig <- unlist(lapply(asplit, function(x)x[1]))
    names(ambig) <- unlist(lapply(asplit, function(x)x[2]))
    ambig <- ambig[ambig != "GERM"]
    nambig <- states[!states %in% names(ambig)]
    names(nambig) <- nambig
    ambig <- c(ambig, states)
  }else{
    ambig <- states
  }
  return(ambig)
}

#' Make a parsimony model file
#' 
#' \code{makeModelFile} Filler
#' @param    file          model file name to write.
#' @param    states        vector of states to include in model.
#' @param    constraints   constraints to add to model.
#' @param    exceptions    vector of comma-separated states that are 
#'                         exceptions to constraints
#'
#'
#' @return   Name of model file
#'
#' @details
#' Currently the only option for \code{constraints} is "irrev", which
#' forbids switches moving from left to right in the \code{states} vector.
#'  
#' @seealso \link{readModelFile}, \link{getTrees}, \link{findSwitches}
#'
#' @export
makeModelFile <- function(file, states, constraints=NULL, exceptions=NULL){
  write("#BEGIN", file=file)
  write(length(states), file=file, append=TRUE)
  write("", file=file, append=TRUE)
  write("#STATES", file=file, append=TRUE)
  for(a in states){
    write(paste(a), file=file, append=TRUE)
  }
  write("", file=file, append=TRUE)
  write("#CONSTRAINTS", file=file, append=TRUE)
  if(!is.null(constraints)){
    if(constraints=="irrev"){
      for(i in 1:(length(states)-1)){
        for(j in (i+1):length(states)){
          if(paste0(states[j],",",states[i]) %in% exceptions){
            print(paste("Excepting",paste0(states[j],",",states[i]),
                        "from constraints"))
            next
          }
          write(paste(states[j], states[i], "1000"), 
                file=file, append=TRUE)
        }
      }
    }else{
      for(constraint in constraints){
        cons <- strsplit(constraint,split=",")[[1]]
        write(paste(cons[1], cons[2], "1000"), 
              file=file, append=TRUE)
      }
    }
  }
  write("", file=file, append=TRUE)
  write("", file=file, append=TRUE)
  write("#AMBIGUOUS", file=file, append=TRUE)
  for(a in states){
    write(paste("GERM", a), file=file, append=TRUE)
  }
  return(file)
}

# Read in a switches file from IgPhyML
# 
# \code{readSwitches} Filler
# @param    file          IgPhyML output file.
#
# @return   A named vector containing the states of the model
readSwitches <- function(file){
  t <- read.table(file,sep="\t",stringsAsFactors=FALSE)
  names(t) <- c("REP","FROM","TO","SWITCHES")
  t
}

# Make bootstrap replicate of clonal alignment
# 
# \code{lones} Filler
# @param    clone     \code{airrClone} object
# @param    reps      Number of bootstrap replicates
# @param    partition If "locus" Bootstrap heavy/lights separately
#
# @return   A list of \code{airrClone} objects with 
# bootstrapped sequences
bootstrapClones  <- function(clone, reps=100, partition="locus", by_codon = TRUE){
  if(clone@phylo_seq == "hlsequence"){
    sarray <- strsplit(clone@data$hlsequence,split="")
    garray <- strsplit(clone@hlgermline,split="")[[1]]
    index <- 1:stats::median(nchar(clone@data$hlsequence))
  }else if(clone@phylo_seq == "sequence"){
    sarray <- strsplit(clone@data$sequence,split="")
    garray <- strsplit(clone@germline,split="")[[1]]
    index <- 1:stats::median(nchar(clone@data$sequence))
  }else if(clone@phylo_seq == "lsequence"){
    sarray <- strsplit(clone@data$lsequence,split="")
    garray <- strsplit(clone@lgermline,split="")[[1]]
    index <- 1:stats::median(nchar(clone@data$lsequence))
  }else{
    stop(paste("phyloseq option",clone@phylo_seq,"not recognized"))
  }
  bootstraps <- list()
  for(i in 1:reps){
    clone_copy <- clone
    if(by_codon == TRUE){
      if(partition == "locus"){
        sindex <- unlist(lapply(unique(clone@locus), function(x){
          locus = which(clone@locus == x)
          seq = seq(1,length(locus),by=3)
          codons=lapply(seq, function(x)c(locus[x],locus[x+1],locus[x+2]))
          codons = sample(codons, replace=TRUE)
          return(codons)}))
        
      }else{
        seq = seq(1,length(index),by=3)
        codons=lapply(seq, function(x)c(index[x],index[x+1],index[x+2]))
        codons = sample(codons, replace=TRUE)
        sindex <- unlist(codons)
      }
    }else{
      if(partition == "locus"){
        sindex <- unlist(lapply(unique(clone@locus), function(x){
          locus = which(clone@locus == x)
          seq = seq(1,length(locus),by=1)
          codons=lapply(seq, function(x)c(locus[x]))
          codons = sample(codons, replace=TRUE)
          return(codons)}))
      }else{
        sindex <- sample(index,length(index),replace=TRUE)
      }
    }
    if(clone@phylo_seq == "sequence"){

      clone_copy@data$sequence <- unlist(lapply(sarray,
                                                function(x)paste(x[sindex],collapse="")))
      clone_copy@germline <- paste(garray[sindex],collapse="")
      if(by_codon == TRUE){
        ## add check step where you translate it and search for \\*
        check <- grepl("\\*", alakazam::translateDNA(clone_copy@data$sequence))
        if(TRUE %in% check){
          trues <- c()
          for(i in 1:length(check)){
            if(check[i] == TRUE){
              trues <- append(trues, clone_copy@data$sequence_id[i])
            }
          }
          stop(paste0("Stop codon(s) have been found in sequence(s) ", trues, " of clone ", clone_copy@clone))
        }
      }
    }else if(clone@phylo_seq == "hlsequence"){
      clone_copy@data$hlsequence <- unlist(lapply(sarray,
                                                  function(x)paste(x[sindex],collapse="")))
      clone_copy@hlgermline <- paste(garray[sindex],collapse="")
      if(by_codon == TRUE){
        ## add check step where you translate it and search for \\*
        check <- grepl("\\*", alakazam::translateDNA(clone_copy@data$hlsequence))
        if(TRUE %in% check){
          trues <- c()
          for(i in 1:length(check)){
            if(check[i] == TRUE){
              trues <- append(trues, clone_copy@data$sequence_id[i])
            }
          }
          stop(paste0("Stop codon(s) have been found in sequence(s) ", trues, " of clone ", clone_copy@clone))
        }
      }
    }else if(clone@phylo_seq == "lsequence"){
      clone_copy@data$lsequence <- unlist(lapply(sarray,
                                                 function(x)paste(x[sindex],collapse="")))
      clone_copy@lgermline <- paste(garray[sindex],collapse="")
      if(by_codon == TRUE){
        ## add check step where you translate it and search for \\*
        check <- grepl("\\*", alakazam::translateDNA(clone_copy@data$lsequence))
        if(TRUE %in% check){
          trues <- c()
          for(i in 1:length(check)){
            if(check[i] == TRUE){
              trues <- append(trues, clone_copy@data$sequence_id[i])
            }
          }
          stop(paste0("Stop codon(s) have been found in sequence(s) ", trues, " of clone ", clone_copy@clone))
        }
      }
    }else{
      stop(paste("phylo_seq",clone@phylo_seq,"not recognized"))
    }
    bootstraps[[i]] <- clone_copy
  }
  bootstraps
}

#' Do IgPhyML maximum parsimony reconstruction
#' 
#' \code{reconIgPhyML} IgPhyML parsimony reconstruction function
#' @param    file          IgPhyML lineage file (see writeLineageFile)
#' @param    modelfile      File specifying parsimony model
#' @param    id            id for IgPhyML run
#' @param    igphyml       location of igphyml executable
#' @param    mode          return trees or count switches? (switches or trees)
#' @param    type          get observed switches or permuted switches?
#' @param    nproc         cores to use for parallelization
#' @param    quiet         amount of rubbish to print
#' @param    rm_files      remove temporary files?
#' @param    rm_dir        remove temporary directory?
#' @param    states        states in parsimony model
#' @param    palette       deprecated
#' @param    resolve       level of polytomy resolution. 0=none, 
#'                         1=maximum parsimony, 2=maximum ambiguity
#' @param    rseed         random number seed if desired
#' @param    force_resolve continue even if polytomy resolution fails?
#' @param    ...           additional arguments
#'
#' @return   Either a tibble of switch counts or a list
#'           of trees with internal nodes predicted by parsimony.
#' @export
reconIgPhyML <- function(file, modelfile, id, 
                         igphyml="igphyml", mode="switches", type="recon",
                         nproc=1, quiet=0, rm_files=FALSE, rm_dir=NULL, 
                         states=NULL, palette=NULL, resolve=2, rseed=NULL,
                         force_resolve=FALSE, ...){
  
  #args <- list(...)
  igphyml <- path.expand(igphyml)
  if(file.access(igphyml, mode=1) == -1) {
    stop("The file ", igphyml, " cannot be executed.")
  }
  if(!file.exists(file)) {
    stop("The repertoire file ", file, " cannot be found.")
  }
  if(!file.exists(modelfile)) {
    stop("The model file ", modelfile, " cannot be found.")
  }
  if(!mode %in% c("switches","trees")){
    stop(paste("mode must be either switches or trees"))
  }
  if(!type %in% c("recon","permute","permuteAll")){
    stop(paste("type must be either recon, permute, or permuteAll"))
  }
  
  if(quiet > 0){
    print(paste("Resolve:",resolve,"Force resolve?",force_resolve))
  }
  recon <- paste0(file,"_igphyml_parstats_",type,".txt")
  logfile <- paste0(file,".log")
  log <- paste(">>",logfile)
  permute <- ""
  force_resolve_option <- ""
  if(type == "permute"){
    permute <- "--permute"
  }
  if(type == "permuteAll"){
    permute <- "--permuteAll"
  }
  if(force_resolve){
    force_resolve_option <- "--force_resolve"
  }
  if(is.null(rseed)){
    rseed <- ""
  }else{
    rseed <- paste("--r_seed",rseed)
  }
  command <- paste("--repfile",file,
                   "--recon",modelfile,"--threads",nproc,"--polyresolve",resolve,
                   "-m HLP -o n --motifs WRC_2:0 --hotness 0 --run_id",type,permute,
                   force_resolve_option,rseed,log)
  params <- list(igphyml,command,stdout=TRUE,stderr=TRUE)
  if(quiet > 2){
    print(paste(params,collapse=" "))
  }
  status <- tryCatch(do.call(base::system2, params), error=function(e){
    print(paste("igphyml error, trying again: ",e));
    cat(paste(readLines(logfile),"\n"))
    return(e)
  }, warning=function(w){
    print(paste("igphyml warnings, trying again: ",w));
    cat(paste(readLines(logfile),"\n"))
    return(w)
  })
  if(length(status) != 0){
    status <- tryCatch(do.call(base::system2, params), error=function(e){
      print(paste("igphyml error, again! quitting: ",e));
      cat(paste(readLines(logfile),"\n"))
      stop()
    }, warning=function(w){
      print(paste("igphyml warnings, again! quitting: ",w));
      cat(paste(readLines(logfile),"\n"))
      stop()
    })
  }
  
  if(quiet > 2){
    cat(paste(readLines(logfile),"\n"))
  }
  if(mode == "switches"){
    recons <- readSwitches(recon)
    recons$CLONE <- id
    recons$TYPE <- toupper(type)
    results <- dplyr::as_tibble(recons)
  }else{
    if(is.null(states)){
      states <- readModelFile(modelfile)
    }
    if(!is.null(palette)){
      warning("palette option is deprecated in reconIgPhyML")
      palette <- NULL
    }
    results <- readLineages(file,states,palette,"recon",quiet)
    results <- lapply(results,function(x){
      x$pars_recon="igphyml";
      x})
    results <- lapply(results,function(x){
      x$tip.label;
      x})
  }
  if(rm_files){
    lines <- readLines(file)
    for(i in 2:length(lines)){
      temp <- strsplit(lines[i],split="\t")[[1]]
      unlink(paste0(temp[1],"*"))
      unlink(paste0(temp[2],"*"))
    }
    unlink(paste0(file,"*"))
  }
  if(!is.null(rm_dir)){
    if(quiet > 1){
      print(paste("rming dir",rm_dir))
    }
    unlink(rm_dir,recursive=TRUE)
  }
  return(results)
}

#' Read in all trees from a lineages file
#' 
#' @param    file    IgPhyML lineage file
#' @param    states  states in parsimony model
#' @param    palette deprecated
#' @param    run_id  id used for IgPhyML run
#' @param    quiet   avoid printing rubbish on screen?
#' @param    append  string appended to fasta files
#' @param    format  format of input file with trees
#' @param    type    Read in parsimony reconstructions or ancestral sequence
#'                   reconstructions? "jointpars" reads in parsimony states, 
#'                   others read in sequences in internal nodes
#'
#' @return   A list of phylo objects from \code{file}.
#' @export
readLineages <- function(file, states=NULL, palette=NULL,
                         run_id="", quiet=TRUE, append=NULL, format="nexus", 
                         type="jointpars"){
  if(!is.null(palette)){
    warning("palette is deprecated in readLineages")
    palette <- NULL
  }
  trees <- list()
  t <- readLines(file)
  if(length(t) == 0){
    return(list())
  }
  for(i in 2:length(t)){
    fasta <- strsplit(t[i],split="\t")[[1]][1]
    if(quiet > 0){print(fasta)}
    if(is.null(append)){
      if(run_id == ""){
        tf <- suppressWarnings(phylotate::read_annotated(
          paste0(fasta,"_igphyml_",type,".nex"),
          format=format))
      }else{
        tf <- suppressWarnings(phylotate::read_annotated(
          paste0(fasta,"_igphyml_",type,"_",run_id,".nex"),
          format=format))
      }
    }else{
      tf <- suppressWarnings(read_annotated(paste0(fasta,append),
                                            format=format))
    }
    if(!is.null(states)){
      tf <- condenseTrees(tf,states,palette)
    }
    germ <- tf$tip.label[grep("_GERM",tf$tip.label)]
    tf$name <- gsub("_GERM$","",germ)
    tf$tip.label[which(tf$tip.label == germ)] <- "Germline"
    nnodes <- length(unique(c(tf$edge[,1],tf$edge[,2])))
    tf$nodes <- rep(list(sequence=NULL),times=nnodes)
    if(type=="jointpars"){
      tf <- rerootTree(tf,"Germline",verbose=0)
    }else{
      tf$nodes <- lapply(1:length(tf$nodes),function(x){
        tf$nodes[[x]]$sequence <- tf$node.comment[x]
        tf$nodes[[x]]
      })
    }
    trees[[i-1]] <- tf
  }
  return(trees)
}


#' Write lineage file for IgPhyML use
#' 
#' @param    data      list of \code{airrClone} objects
#' @param    trees     list of \code{phylo} objects corresponding to \code{data}
#' @param    dir       directory to write file
#' @param    id        id used for IgPhyML run
#' @param    rep       bootstrap replicate
#' @param    trait     string appended to sequence id in fasta files
#' @param    partition how to partition omegas
#' @param    heavy     name of heavy chain locus
#' @param    empty     output uninformative sequences?
#' @return   Name of created lineage file.
#' @export
writeLineageFile <- function(data, trees=NULL, dir=".", id="N", rep=NULL, 
                             trait=NULL, empty=TRUE, partition="single", heavy="IGH"){
  
  file <- file.path(dir,paste0(id,"_lineages_pars.tsv"))
  if(!is.null(rep)){
    file <- file.path(dir,paste0(id,"_lineages_",rep,"_pars.tsv"))
  }
  outdir <- file.path(dir,paste0(id,"_recon_",rep))
  dir.create(dir,showWarnings=FALSE)
  dir.create(outdir,showWarnings=FALSE)
  
  dnames <- unlist(lapply(data,function(x)x@clone))
  tnames <- unlist(lapply(trees,function(x)x$name))
  if(sum(tnames != dnames) != 0){
    trees <- trees[order(match(tnames,dnames))]
  }
  write(length(data),file=file)
  for(i in 1:length(data)){
    tree <- trees[[i]]
    fastafile <- file.path(outdir,paste0(data[[i]]@clone,".fasta"))
    treefile <- file.path(outdir,paste0(data[[i]]@clone,".tree"))
    germid <- paste0(data[[i]]@clone,"_GERM")
    cloneToFasta(data[[i]],fastafile,germid,trait,empty=empty)
    if(data[[i]]@phylo_seq == "sequence"){
      g <- data[[i]]@germline
    }else if(data[[i]]@phylo_seq == "lsequence"){
      g <- data[[i]]@lgermline
    }else if(data[[i]]@phylo_seq == "hlsequence"){
      g <- data[[i]]@hlgermline
    }else{
      stop(paste("phylo_seq not recognized",c@clone))
    }
    
    if(!partition %in% c("single","hl")){
      acceptable <- c("fwr1","fwr2","fwr3","fwr4","cdr1","cdr2","cdr3")
      unacceptable <- unlist(lapply(data, function(x)sum(!x@region %in% acceptable) > 0))
      exclude_clones <- unlist(lapply(data[unacceptable], function(x)x@clone))
      if(length(exclude_clones) > 0){
        stop(paste("non-standard regions found in these clones,",
                   "either remove or set partition='single':",paste(exclude_clones, collapse=","),
                   "\nAllowable regions:",paste(acceptable,collapse=",")))
      }
    }
    
    if(partition == "cf"){ #make file specifying sequence regions
      nomega <- 2
      regions <- data[[i]]@region
      if(dplyr::n_distinct(regions) == 1){
        warning(paste("Only one region found in clone",data[[i]]@clone))
      }
      cdrs <- rep(0,length(regions))
      cdrs[regions == "fwr1"] <- 13
      cdrs[regions == "cdr1"] <- 30
      cdrs[regions == "fwr2"] <- 45
      cdrs[regions == "cdr2"] <- 60
      cdrs[regions == "fwr3"] <- 80
      cdrs[regions == "cdr3"] <- 108
      cdrs[regions == "fwr4"] <- 120
    }else if(partition == "hl"){
      nomega <- 2
      chains <- data[[i]]@locus
      if(dplyr::n_distinct(chains) == 1){
        warning(paste("Only one chain found in clone",data[[i]]@clone))
      }
      cdrs <- rep(0,length(chains))
      cdrs[chains == heavy] <- 13
      cdrs[chains != heavy] <- 30
    }else if(partition == "hlc"){
      nomega <- 3
      chains <- data[[i]]@locus
      regions <- data[[i]]@region
      if(dplyr::n_distinct(regions) == 1){
        warning(paste("Only one region found in clone",data[[i]]@clone))
      }
      if(dplyr::n_distinct(chains) == 1){
        warning(paste("Only one chain found in clone",data[[i]]@clone))
      }
      cdrs <- rep(0,length(chains))
      cdrs[grepl("fwr", regions)] <- 13
      cdrs[chains == heavy & grepl("cdr", regions)] <- 30  #heavy cdr
      cdrs[chains != heavy & grepl("cdr", regions)] <- 200 #light cdr
    }else if(partition == "hlf"){
      nomega <- 3
      chains <- data[[i]]@locus
      regions <- data[[i]]@region
      if(dplyr::n_distinct(regions) == 1){
        warning(paste("Only one region found in clone",data[[i]]@clone))
      }
      if(dplyr::n_distinct(chains) == 1){
        warning(paste("Only one chain found in clone",data[[i]]@clone))
      }
      cdrs <- rep(0,length(chains))
      cdrs[chains == heavy & grepl("fwr", regions)] <- 13 #heavy fwr
      cdrs[grepl("cdr", regions)] <- 30
      cdrs[chains != heavy & grepl("fwr", regions)] <- 200 #light fwr
    }else if(partition == "hlcf"){
      nomega <- 4
      chains <- data[[i]]@locus
      regions <- data[[i]]@region
      if(dplyr::n_distinct(regions) == 1){
        warning(paste("Only one region found in clone",data[[i]]@clone))
      }
      if(dplyr::n_distinct(chains) == 1){
        warning(paste("Only one chain found in clone",data[[i]]@clone))
      }
      cdrs <- rep(0,length(chains))
      cdrs[chains == heavy & grepl("fwr", regions)] <- 13
      cdrs[chains == heavy & grepl("cdr", regions)] <- 30
      cdrs[chains != heavy & grepl("fwr", regions)] <- 200 #light fwr
      cdrs[chains != heavy & grepl("cdr", regions)] <- 300 #light cdr
    }else if(partition != "single"){
      stop(paste("Partition",partition,"not recognized"))
    }
    
    if(partition != "single"){
      partfile <- file.path(outdir,paste0(data[[i]]@clone,".part.txt"))
      write(paste(nomega,nchar(g)/3), file=partfile)
      write("FWR:IMGT\nCDR:IMGT", file=partfile, append=TRUE)
      if(partition == "hlf" || partition == "hlcf"){
        write("FWRL:IMGT", file=partfile, append=TRUE)
      }
      if(partition == "hlc" || partition == "hlcf"){
        write("CDRL:IMGT", file=partfile, append=TRUE)
      }
      write(paste(data[[i]]@v_gene,data[[i]]@j_gene,sep="\n"), file=partfile, append=TRUE)
      write(paste(cdrs[1:length(cdrs) %% 3 == 0],collapse=","), file=partfile, append=TRUE)
    }else{
      partfile <- "N"
    }
    
    if(!is.null(trees)){
      if("node.label" %in% names(tree)){
        tree$node.label <- NULL
      }
      tree$tip.label[which(tree$tip.label == "Germline")] <- germid
      tree <- ape::multi2di(tree)
      ape::write.tree(tree,file=treefile)
      write(paste(fastafile,treefile,germid,partfile,sep="\t"), file=file,
            append=TRUE)
    }else{
      write(paste(fastafile,"N",germid,partfile,sep="\t"), file=file,
            append=TRUE)
    }
  }
  return(file)
}

#' Wrapper for alakazam::buildPhylipLineage
#' 
#' @param    clone      \code{airrClone} object
#' @param    exec       dnapars or dnaml executable
#' @param    temp_path  path to temporary directory
#' @param    verbose    amount of rubbish to print
#' @param    rm_temp    remove temporary files?
#' @param    seq        sequece column in \code{airrClone} object
#' @param    tree       fixed tree topology if desired (currently does nothing
#'                      if specified)
#' @param    onetree    Only sample one tree if multiple found.
#'
#' @return  \code{phylo} object created by dnapars or dnaml with nodes attribute
#'          containing reconstructed sequences.
#' @export
buildPhylo <- function(clone, exec, temp_path=NULL, verbose=0,
                       rm_temp=TRUE, seq="sequence", tree=NULL, onetree=TRUE){
  
  if(grepl("dnaml$",exec) | grepl("dnaml\\.exe$",exec)){
    method <- "dnaml"
  }else if(grepl("dnapars$",exec) | grepl("dnapars\\.exe$",exec)){
    method <- "dnapars"
  }else{
    stop("executable not recognized! Must end with dnapars or dnaml")
  }
  if(seq != "sequence"){
    clone@data$sequence <- clone@data[[seq]]
    if(seq == "hlsequence"){
      clone@germline <- clone@hlgermline
    }else if(seq == "lsequence"){
      clone@germline <- clone@lgermline
    }
  }
  if(nrow(clone@data) < 2){
    stop("Clone ",paste0(clone@clone," has only one sequence, skipping"))
  }
    if(is.null(tree)){
        tree <- tryCatch(
            alakazam::buildPhylipLineage(clone, exec, rm_temp=rm_temp,
                branch_length="distance", verbose=verbose>0,
                temp_path=temp_path,onetree=onetree),
            error=function(e)e)
        if(is.null(tree)){
            stop(paste0("buildPhylipLineage failed for clone ",clone@clone))
        }
        if(inherits(tree, "error")){
            stop(tree)
        }
    tree <- alakazam::graphToPhylo(tree)
    tree <- rerootTree(tree, germline="Germline",verbose=0)
    tree <- ape::ladderize(tree,right=FALSE)
    tree$name <- clone@clone
    tree$tree_method <- paste0("phylip::",method)
    tree$edge_type <- "genetic_distance"
    tree$seq <- seq
  }
  tree
}

#' Wrapper for phangorn::pratchet
#' 
#' @param    clone           \code{airrClone} object
#' @param    seq             sequece column in \code{airrClone} object
#' @param    asr             return sequence or probability matrix?
#' @param    asr_thresh      threshold for including a nucleotide as an alternative
#' @param    tree            fixed tree topology if desired.
#' @param    asr_type        MPR or ACCTRAN
#' @param    verbose         amount of rubbish to print
#' @param    resolve_random  randomly resolve polytomies?
#' @param    data_type       Are sequences DNA or AA?
#' @return  \code{phylo} object created by phangorn::pratchet with nodes
#'          attribute containing reconstructed sequences.
#' @export
buildPratchet <- function(clone, seq="sequence", asr="seq", asr_thresh=0.05, 
                          tree=NULL, asr_type="MPR", verbose=0, resolve_random=TRUE,
                          data_type="DNA"){
  seqs <- clone@data[[seq]]
  names <- clone@data$sequence_id

  if(verbose > 0){
    print(clone@clone)
  }
  if(length(seqs) < 2){
    stop(paste0(clone@clone," has only one sequence, skipping"))
  }
  if(seq == "hlsequence"){
    germline <- clone@hlgermline
  }else if(seq == "lsequence"){
    germline <- clone@lgermline
  }else{
    germline <- clone@germline
  }
  
  seqs <- base::append(seqs,germline)
  names <- c(names,"Germline")
  
  seqs <- strsplit(seqs,split="")
  names(seqs) <- names

  lengths = unlist(lapply(seqs,function(x)length(x)))
  if(any(lengths != lengths[1])){
    stop(paste0("Sequence and/or germline lengths of clone ",
                clone@clone," are not equal"))
  }
  if(data_type=="DNA"){
    data <- phangorn::phyDat(ape::as.DNAbin(t(as.matrix(dplyr::bind_rows(seqs)))))
  }else{
    data <- phangorn::phyDat(ape::as.AAbin(t(as.matrix(dplyr::bind_rows(seqs)))),
                             type="AA")
  }
  if(is.null(tree)){
    tree <- tryCatch(phangorn::pratchet(data,trace=FALSE),warning=function(w)w)
    tree <- phangorn::acctran(ape::multi2di(tree,random=resolve_random),data)
    tree <- ape::unroot(tree)
    tree$edge.length <- tree$edge.length/nchar(germline)
    tree$tree_method <- "phangorn::prachet"
    tree$edge_type <- "genetic_distance"
    nnodes <- length(unique(c(tree$edge[,1],tree$edge[,2])))
    tree$nodes <- rep(list(sequence=NULL),times=nnodes)
    tree$node.label <- NULL
  }else if(is.null(tree$nodes)){
    nnodes <- length(unique(c(tree$edge[,1],tree$edge[,2])))
    tree$nodes <- rep(list(sequence=NULL),times=nnodes)
  }
  tree$name <- clone@clone
  tree$seq <- seq
  
  if(asr != "none" && data_type=="DNA"){
    seqs_pars <- phangorn::ancestral.pars(tree, data, 
                                          type=asr_type, cost=NULL, return="prob")
    
    # CGJ 5/23/24 for the phangorn update (3.0.0) and so it passes win builder devel
    if("prob" %in% names(seqs_pars)){
      ASR <- list()
      tip_data <- data.frame(seqs_pars$data)
      for(i in 1:ncol(tip_data)){
        if(asr == "seq"){
          sub_seq <- paste0(toupper(tip_data[,i]))
          seq_ar <- unlist(lapply(1:length(sub_seq), function(x){
            site <- sub_seq[x]
            site <- alakazam::DNA_IUPAC[[paste(sort(site), collapse = "")]]
            if(length(site) == 0){
              site <- "N"
            }
            site}))
          ASR[[as.character(i)]] <- paste(seq_ar,collapse="")
        } else{
          ASR[[as.character(i)]] <- paste0(toupper(tip_data[,i]), collapse = "")
        }
      }
      nASR <- length(ASR)
      internal_seqs <- seqs_pars$prob
      for(i in 1:length(unique(internal_seqs$Node))){
        sub <- internal_seqs[internal_seqs$Node == unique(internal_seqs$Node)[i],]
        if(asr == "seq"){
          sub_seq <- paste0(toupper(sub$State))
          seq_ar <- unlist(lapply(1:length(sub_seq), function(x){
            site <- sub_seq[x]
            site <- alakazam::DNA_IUPAC[[paste(sort(site), collapse = "")]]
            if(length(site) == 0){
              site <- "N"
            }
            site}))
          ASR[[as.character(nASR + i)]] <- paste(seq_ar,collapse="")
        } else{
          ASR[[as.character(nASR + i)]] <- paste0(toupper(sub$State), collapse = "")
        }
      }
    } else{ # if they are using an older version of phangorn keep it the same
      ASR <- list()
      for(i in 1:max(tree$edge)){
        patterns <- t(subset(seqs_pars, i)[[1]])
        pat <- patterns[,attr(seqs_pars,"index")]
        if(asr == "seq"){
          thresh <- pat > asr_thresh
          acgt <- c("A","C","G","T")
          seq_ar <- unlist(lapply(1:ncol(pat),function(x){
            site <- acgt[thresh[,x]]
            site <- alakazam::DNA_IUPAC[[paste(sort(site),collapse="")]]
            if(length(site) == 0){
              site <- "N"
            }
            site}))
          ASR[[as.character(i)]] <- paste(seq_ar,collapse="")
        }else{
          ASR[[as.character(i)]] <- pat
        }
      }
    }
    

    tree$nodes <- lapply(1:length(tree$nodes),function(x){
      tree$nodes[[x]]$sequence <- ASR[[x]]
      tree$nodes[[x]]
    })
  }
  opars <- phangorn::parsimony(ape::di2multi(tree),data)
  tree <- rerootTree(tree,"Germline",verbose=0)
  npars <- phangorn::parsimony(ape::di2multi(tree),data)
  if(npars != opars){
    stop(paste("Error in rerooting tree",tree$name,
               "parsimony score not consistent:",opars,npars))
  }
  return(tree)
}


#' Wrapper for phangorn::optim.pml
#' 
#' @param    clone      \code{airrClone} object
#' @param    seq        sequece column in \code{airrClone} object
#' @param    sub_model  substitution model to use
#' @param    gamma      gamma site rate variation?
#' @param    asr        return sequence or probability matrix?
#' @param    asr_thresh threshold for including a nucleotide as an alternative
#' @param    tree       fixed tree topology if desired.
#' @param    data_type  Are sequences DNA or AA?
#' @param    verbose    Print error messages as they happen?
#' @param    optNni     Optimize tree topology
#' @param    optQ       Optimize Q matrix
#' @param    resolve_random  randomly resolve polytomies?
#' @param    quiet           amount of rubbish to print to console
#' @param    rep             current bootstrap replicate (experimental)
#'
#' @return  \code{phylo} object created by phangorn::optim.pml with nodes
#'          attribute containing reconstructed sequences.
#' @export
buildPML <- function(clone, seq="sequence", sub_model="GTR", gamma=FALSE, asr="seq", 
                     asr_thresh=0.05, tree=NULL, data_type="DNA", optNni=TRUE, optQ=TRUE, 
                     verbose=FALSE, resolve_random=TRUE, quiet=0, rep=NULL){
  seqs <- clone@data[[seq]]
  names <- clone@data$sequence_id
  if(length(seqs) < 2){
    stop(paste0(clone@clone," has only one sequence, skipping"))
  }
  if(seq == "hlsequence"){
    germline <- clone@hlgermline
  }else if(seq == "lsequence"){
    germline <- clone@lgermline
  }else{
    germline <- clone@germline
  }
  seqs <- base::append(seqs,germline)
  names <- c(names,"Germline")
  seqs <- strsplit(seqs,split="")
  names(seqs) <- names
  if(data_type=="DNA"){
    data <- phangorn::phyDat(ape::as.DNAbin(t(as.matrix(dplyr::bind_rows(seqs)))))
  }else{
    data <- phangorn::phyDat(ape::as.AAbin(t(as.matrix(dplyr::bind_rows( seqs)))),
                             type="AA")
    if(sub_model == "GTR"){
      warning("GTR model shouldn't be used for AA.")
    }
  }
  if(is.null(tree)){
    dm  <- phangorn::dist.ml(data)
    treeNJ  <- ape::multi2di(phangorn::NJ(dm), random=resolve_random)
    treeNJ$edge.length[treeNJ$edge.length < 0] <- 0 #change negative edge lengths to zero
    pml <- phangorn::pml(ape::unroot(treeNJ),data=data)
    fit <- tryCatch(phangorn::optim.pml(pml, model=sub_model, optNni=optNni, optQ=optQ,
                                        optGamma=gamma, rearrangement="NNI",control=phangorn::pml.control(epsilon=1e-08,
                                        maxit=10, trace=0)), error=function(e)e)
    if("error" %in% class(fit)){
      if(verbose){
        print(fit)
      }
      return(fit)
    }
    #tree <- ape::unroot(ape::multi2di(fit$tree)) #CGJ 4/5/23
    # this assumes we can change tree object without affecting ASR
    fit$tree <- ape::unroot(ape::multi2di(fit$tree))
    tree <- fit$tree
    # test if the tree is binary 
    if(!ape::is.binary(tree)){
      stop(paste("Tree may not be full resolved at", clone@clone))
    }
    tree$tree_method <- paste("phangorn::optim.pml::",sub_model)
    tree$edge_type <- "genetic_distance"
    nnodes <- length(unique(c(tree$edge[,1],tree$edge[,2])))
    tree$nodes <- rep(list(sequence=NULL),times=nnodes)
  }
  tree$name <- clone@clone
  tree$seq <- seq

  if(asr != "none" && data_type=="DNA"){
    seqs_ml <- phangorn::ancestral.pml(fit,
                                       type="marginal",return="prob")
    # CGJ 5/23/24 for the phangorn update (3.0.0) and so it passes win builder devel
    if("prob" %in% names(seqs_ml)){
      ASR <- list()
      tip_data <- data.frame(seqs_ml$data)
      for(i in 1:ncol(tip_data)){
        if(asr == "seq"){
          sub_seq <- paste0(toupper(tip_data[,i]))
          seq_ar <- unlist(lapply(1:length(sub_seq), function(x){
            site <- sub_seq[x]
            site <- alakazam::DNA_IUPAC[[paste(sort(site), collapse = "")]]
            if(length(site) == 0){
              site <- "N"
            }
            site}))
          ASR[[as.character(i)]] <- paste(seq_ar,collapse="")
        } else{
          ASR[[as.character(i)]] <- paste0(toupper(tip_data[,i]), collapse = "")
        }
      }
      nASR <- length(ASR)
      internal_seqs <- seqs_ml$prob
      for(i in 1:length(unique(internal_seqs$Node))){
        sub <- internal_seqs[internal_seqs$Node == unique(internal_seqs$Node)[i],]
        if(asr == "seq"){
          sub_seq <- paste0(toupper(sub$State))
          seq_ar <- unlist(lapply(1:length(sub_seq), function(x){
            site <- sub_seq[x]
            site <- alakazam::DNA_IUPAC[[paste(sort(site), collapse = "")]]
            if(length(site) == 0){
              site <- "N"
            }
            site}))
          ASR[[as.character(nASR + i)]] <- paste(seq_ar,collapse="")
        } else{
          ASR[[as.character(nASR + i)]] <- paste0(toupper(sub$State), collapse = "")
        }
      }
    } else{ # if they are using an older version of phangorn keep it the same
      ASR <- list()
      for(i in 1:max(tree$edge)){
        patterns <- t(subset(seqs_ml, i)[[1]])
        pat <- patterns[,attr(seqs_ml,"index")]
        if(asr == "seq"){
          thresh <- pat > asr_thresh
          acgt <- c("A","C","G","T")
          seq_ar <- unlist(lapply(1:ncol(pat),function(x){
            site <- acgt[thresh[,x]]
            site <- alakazam::DNA_IUPAC[[paste(sort(site),collapse="")]]
            if(length(site) == 0){
              site <- "N"
            }
            site}))
          ASR[[as.character(i)]] <- paste(seq_ar,collapse="")
        }else{
          ASR[[as.character(i)]] <- pat
        }
      }
    }
    tree$nodes <- lapply(1:length(tree$nodes),function(x){
      tree$nodes[[x]]$sequence <- ASR[[x]]
      tree$nodes[[x]]
    })
  }
  tree <- rerootTree(tree,"Germline",verbose=0)
  tree$parameters <- fit
  return(tree)
}




#' Wrapper to build IgPhyML trees and infer intermediate nodes
#' 
#' @param    clone      list of \code{airrClone} objects
#' @param    igphyml    igphyml executable
#' @param    trees      list of tree topologies if desired
#' @param    nproc      number of cores for parallelization
#' @param    temp_path  path to temporary directory
#' @param    id         IgPhyML run id
#' @param    rseed      random number seed if desired
#' @param    quiet      amount of rubbish to print
#' @param    rm_files   remove temporary files?
#' @param    rm_dir     remove temporary directory?
#' @param    partition  How to partition omegas along sequences (see details)
#' @param    omega      omega parameters to estimate (see IgPhyML docs)
#' @param    optimize   optimize HLP rates (r), lengths (l), topology (t)
#' @param    motifs     motifs to consider (see IgPhyML docs)
#' @param    hotness    hotness parameters to estimate (see IgPhyML docs)
#' @param    rates      comma delimited list showing which omega-defined partitions
#'                      get a separate rate (e.g. omega=e,e rates=0,1).
#' @param    asrc       Intermediate sequence cutoff probability
#' @param    splitfreqs Calculate codon frequencies on each partition separately?
#' @param    ...        Additional arguments (not currently used)
#'
#' @details Partition options in rate order:
#' \itemize{
#'   \item  \code{single}: 1 omega for whole sequence
#'   \item  \code{cf}: 2 omegas, 1 for all CDRs and 1 for all FWRs
#'   \item  \code{hl}: 2 omegas, 1 for heavy and 1 for light chain
#'   \item  \code{hlf}: 3 omegas, 1 for heavy FWR, 1 for all CDRs, and 1 for light FWRs
#'   \item  \code{hlc}: 3 omegas, 1 for all FWRs, 1 for heavy CDRs, and 1 for light CDRs
#'   \item  \code{hlcf}: 4 omegas, 1 for each heavy FWR, 1 for heavy CDR, 1 for light FWR, and 1 for light CDR
#' }
#'
#' @return  \code{phylo} object created by igphyml with nodes attribute
#'          containing reconstructed sequences.
#' @export
buildIgphyml <- function(clone, igphyml, trees=NULL, nproc=1, temp_path=NULL, 
                         id=NULL, rseed=NULL, quiet=0, rm_files=TRUE, rm_dir=NULL, 
                         partition=c("single", "cf", "hl", "hlf", "hlc", "hlcf"),
                         omega=NULL, optimize="lr", motifs="FCH", hotness="e,e,e,e,e,e", 
                         rates=NULL, asrc=0.95, splitfreqs=FALSE, ...){
  
  warning("Dowser igphyml doesn't mask split codons!")
  partition <- match.arg(partition)
  
  valid_o <- c("r","lr","tlr")
  if(!optimize %in% valid_o){
    stop(paste("Invalid optimize specification, must be one of:",valid_o))
  }
  na_regions <- unlist(lapply(clone, function(x)sum(is.na(x@region)) > 0))
  if(sum(na_regions) > 0){
    exclude_clones <- unlist(lapply(clone[na_regions], function(x)x@clone))
    stop(paste("NA regions found in clones",paste(exclude_clones, collapse=","), 
               "remove before continuing"))
  }
  
  if(!is.null(omega)){
    os <- strsplit(omega,split=",")[[1]]
  } else {
    os <- NULL
  }
  file <- writeLineageFile(clone,trees,dir=temp_path,id=id,rep=id,empty=FALSE,
                           partition=partition, ...)
  if(length(os) != 2 && (partition == "cf" | partition == "hl")){
    warning("Omega parameter incompatible with partition, setting to e,e")
    omega = "e,e"
    os <- strsplit(omega,split=",")[[1]]
    if(partition == "hl" && is.null(rates)){
      rates = "0,1"
    }
  }
  if(length(os) != 3 && (partition == "hlc" | partition == "hlf")){
    warning("Omega parameter incompatible with partition, setting to e,e,e")
    omega = "e,e,e"
    os <- strsplit(omega,split=",")[[1]]
  }
  if(length(os) != 4 && (partition == "hlcf")){
    warning("Omega parameter incompatible with partition, setting to e,e,e,e")
    omega = "e,e,e,e"
    os <- strsplit(omega,split=",")[[1]]
    if(is.null(rates)){
      rates = "0,0,1,1"
    }
  }
  if(length(os) != 1 && (partition == "single")){
    warning("Omega parameter incompatible with partition, setting to e")
    omega = "e"
    os <- strsplit(omega,split=",")[[1]]
    #stop("Specified partition model not compatible with multiple omegas or rates")
  }
  if (!is.null(rates)) {
    if (length(os) != length(strsplit(rates, ",")[[1]])) {
      stop("Number of rates needs to equal the number of omegas")
    }
  }
  igphyml <- path.expand(igphyml)
  if(file.access(igphyml, mode=1) == -1) {
    stop("The file ", igphyml, " cannot be executed.")
  }
  logfile <- paste0(file,".log")
  log <- paste(">>",logfile)
  if(is.null(rseed)){
    rseed <- ""
  }else{
    rseed <- paste("--r_seed",rseed)
  }
  gyrep <- paste0(file,"_gyrep")
  if(!is.null(trees)){
    command <- paste("--repfile",file,"--outrep",gyrep,
                     "--threads",nproc,"-o lr -m GY --run_id gy",rseed,log)
  }else{
    command <- paste("--repfile",file,"--outrep",gyrep,
                     "--threads",nproc,"-o tlr -m GY --run_id gy",rseed,log)
  }
  params <- list(igphyml,command,stdout=TRUE,stderr=TRUE)
  if(quiet > 2){
    print(paste(params,collapse=" "))
  }
  status <- tryCatch(do.call(base::system2, params), error=function(e){
    print(paste("igphyml error, trying again: ",e));
    cat(paste(readLines(logfile),"\n"))
    return(e)
  }, warning=function(w){
    print(paste("igphyml warnings, trying again: ",w));
    cat(paste(readLines(logfile),"\n"))
    return(w)
  })
  if(length(status) != 0){
    status <- tryCatch(do.call(base::system2, params), error=function(e){
      print(paste("igphyml error, again! quitting: ",e));
      cat(paste(readLines(logfile),"\n"))
      stop()
    }, warning=function(w){
      print(paste("igphyml warnings, again! quitting: ",w));
      cat(paste(readLines(logfile),"\n"))
      stop()
    })
  }
  if(splitfreqs){
    splitf = "--splitfreqs"
  }else{
    splitf = ""
  }
  if(!is.null(rates)){
    ratestring = paste0("--rates ",rates)
  }else{
    ratestring = ""
  }
  command <- paste("--repfile",gyrep,
                   "--threads",nproc,"--omega",omega,"-o",optimize,"--motifs",motifs,
                   "--hotness",hotness,"-m HLP --run_id hlp --oformat tab --ASRc",asrc,
                   ratestring,splitf,rseed,log)
  params <- list(igphyml,command,stdout=TRUE,stderr=TRUE)
  if(quiet > 2){
    print(paste(params,collapse=" "))
  }
  status <- tryCatch(do.call(base::system2, params), error=function(e){
    print(paste("igphyml error, trying again: ",e));
    cat(paste(readLines(logfile),"\n"))
    return(e)
  }, warning=function(w){
    print(paste("igphyml warnings, trying again: ",w));
    cat(paste(readLines(logfile),"\n"))
    return(w)
  })
  if(length(status) != 0){
    status <- tryCatch(do.call(base::system2, params), error=function(e){
      print(paste("igphyml error, again! quitting: ",e));
      cat(paste(readLines(logfile),"\n") )
      stop()
    }, warning=function(w){
      print(paste("igphyml warnings, again! quitting: ",w));
      cat(paste(readLines(logfile),"\n"))
      stop()
    })
  }
  #trees <- readLineages(file=gyrep,run_id="hlp",type="asr")
  ofile <- file.path(temp_path,paste0(id,"_lineages_",id,
                                      "_pars.tsv_gyrep_igphyml_stats_hlp.tab"))
  results <- alakazam::readIgphyml(ofile,format="phylo",
                                   branches="distance")
  if(partition == "hl"){
    names(results$param) = gsub("fwr","heavy",names(results$param))
    names(results$param) = gsub("cdr","light",names(results$param))
  }else if(partition == "hlf"){
    names(results$param) = gsub("omega_1","omega_heavyfwr",names(results$param))
    names(results$param) = gsub("omega_2","omega_cdr",names(results$param))
    names(results$param) = gsub("omega_3","omega_lightfwr",names(results$param))
  }else if(partition == "hlc"){
    names(results$param) = gsub("omega_1","omega_fwr",names(results$param))
    names(results$param) = gsub("omega_2","omega_heavycdr",names(results$param))
    names(results$param) = gsub("omega_3","omega_lightcdr",names(results$param))
  }else if(partition == "hlcf"){
    names(results$param) = gsub("omega_1","omega_heavyfwr",names(results$param))
    names(results$param) = gsub("omega_2","omega_heavycdr",names(results$param))
    names(results$param) = gsub("omega_3","omega_lightfwr",names(results$param))
    names(results$param) = gsub("omega_4","omega_lightcdr",names(results$param))
  }
  ASR <- readFasta(file.path(temp_path, paste0(id,"_lineages_",id,"_pars_hlp_asr.fasta")))
  trees <- results$trees
  params <- results$param[-1,]
  for(i in 1:nrow(params)){
    clone_id <- params[i,]$clone
    trees[[i]]$name <- clone_id
    trees[[i]]$tip.label[which(trees[[i]]$tip.label 
                               == paste0(clone_id,"_GERM"))] <- "Germline"
    trees[[i]]$tree_method <- paste("igphyml::gy94,hlp19")
    trees[[i]]$edge_type <- "genetic_distance_codon"
    trees[[i]]$seq <- "sequence"
    trees[[i]]$parameters <- c(params[i,])
    # add sequences to internal nodes
    gline_id <- paste0(clone_id,"_GERM")
    ntips <- length(trees[[i]]$tip.label)
    nint  <- dplyr::n_distinct(trees[[i]]$edge[,1])
    nnodes <- ntips + nint
    if(nnodes != length(unique(c(trees[[i]]$edge[,1],
                                 trees[[i]]$edge[,2])))){
      stop(paste("Internal node count error, clone ",clone_id))
    }

    if("rate_light_mle" %in% names(params)){
      #re-scale branches if separate rate estimated
      clone_index <- which(sapply(clone, function(x)x@clone == clone_id))
      if(length(clone_index) != 1){
        stop(paste0("Clone index error: ",clone_id))
      }
      clone_obj <- clone[[clone_index]]
      lrate <- params[i,]$rate_light_mle
      hrate <- params[i,]$rate_heavy_mle
      heavy_sites <- sum(clone_obj@locus == "IGH")/3
      light_sites <- sum(clone_obj@locus != "IGH")/3
      l <- trees[[i]]$edge.length
      trees[[i]]$edge.length <- (l*hrate*heavy_sites + l*lrate*light_sites)/
      (heavy_sites + light_sites)
    }

    trees[[i]]$nodes <- rep(list(sequence=NULL),times=nnodes)
    # blank node should be MRCA, and all labels should be in ASR file
    labs <- trees[[i]]$node.label
    if(sum(labs == "") != 1 || sum(labs == gline_id) > 0){
      stop(paste("Error in reading node labels, clone",clone_id))
    }
    labs[labs == ""] <- gline_id
    if(sum(!labs %in% names(ASR)) != 0){
      stop(paste("Labels not in reconstructed clone",clone_id,":",
                 paste(labs[!labs %in% names(ASR)],collapse=", ")))
    }
    # also add tip sequences to tip nodes
    tipseqs <- clone[[i]]@data[[clone[[i]]@phylo_seq]]
    names(tipseqs) <- clone[[i]]@data$sequence_id
    if(clone[[i]]@phylo_seq == "sequence"){
      tipseqs <- c(tipseqs,"Germline"=clone[[i]]@germline)
    }else if(clone[[i]]@phylo_seq == "lsequence"){
      tipseqs <- c(tipseqs,"Germline"=clone[[i]]@lgermline)
    }else if(clone[[i]]@phylo_seq == "hlsequence"){
      tipseqs <- c(tipseqs,"Germline"=clone[[i]]@hlgermline)
    }else{
      stop(paste("phylo_seq not recognized",clone_id))
    }
    if(sum(!trees[[i]]$tip.label %in% names(tipseqs)) != 0 ||
       sum(!names(tipseqs) %in% trees[[i]]$tip.label) != 0){
      stop(paste("Tip sequences do not match in clone",clone_id))
    }
    trees[[i]]$nodes <- lapply(1:nnodes,function(x){
      if(x <= ntips){ # if a tip, just add sequence
        trees[[i]]$nodes[[x]]$sequence <- tipseqs[trees[[i]]$tip.label[x]]
        trees[[i]]$nodes[[x]]$id <- trees[[i]]$tip.label[x]
      }else{ # if internal node, add reconstructed sequence
        trees[[i]]$nodes[[x]]$sequence <- ASR[[labs[x-ntips]]]
        trees[[i]]$nodes[[x]]$id <- labs[x-ntips]
        if(labs[x-ntips] == gline_id){trees[[i]]$nodes[[x]]$id <- "Germline_Inferred"}
      }
      trees[[i]]$nodes[[x]]
    })
  }

  
  if(rm_files){
    lines <- readLines(file)
    for(i in 2:length(lines)){
      temp <- strsplit(lines[i],split="\t")[[1]]
      unlink(paste0(temp[1],"*"))
      unlink(paste0(temp[2],"*"))
      if(temp[3] != "N"){
        unlink(paste0(temp[3],"*"))
      }
    }
    unlink(paste0(file,"*"))
    unlink(file.path(temp_path,paste0(id,"_lineages_",id,
                                      "_pars_hlp_asr.fasta")))
  }
  if(!is.null(rm_dir)){
    if(quiet > 1){
      print(paste("rming dir",rm_dir))
    }
    unlink(rm_dir,recursive=TRUE)
  }
  return(trees)
}


#' Wrapper to build RAxML-ng trees and infer intermediate nodes
#' 
#' @param    clone      list of \code{airrClone} objects
#' @param    seq        the phylo_seq option does this clone uses. Possible options are "sequence", "hlsequence", or "lsequence"
#' @param    exec       RAxML-ng executable
#' @param    model      The DNA model to be used. GTR is the default.
#' @param    partition  A parameter that determines how branches are reported when partitioning. Options include NULL (default), 
#'                      scaled, unlinked, and linked
#' @param    rseed      The random seed used for the parsimony inferences. This allows you to reproduce your results.
#' @param    name       specifies the name of the output file
#' @param    starting_tree specifies a user starting tree file name and path in Newick format
#' @param    data_type  Specifies what format your data is in, DNA or AA
#' @param    from_getTrees A logical that indicates if the desired starting tree is from getTrees and not a newick file
#' @param    rm_files   remove temporary files?
#' @param    asr        computes the marginal ancestral states of a tree
#' @param    rep        Which repetition of the tree building is currently being run. Mainly for getBootstraps. 
#' @param    dir        Where the output files are to be made. 
#' @param    n_starts   Number of max parsimony starting trees (default is 10 pars + 10 random)
#' @param    ...        Additional arguments (not currently used)
#'
#'
#' @return  \code{phylo} object created by RAxML-ng with nodes attribute
#'          containing reconstructed sequences.
#' @export
buildRAxML <- function(clone, seq = "sequence", exec, model = 'GTR', partition = NULL, 
                       rseed = 28, name = "run", starting_tree = NULL, data_type = "DNA",
                       from_getTrees = FALSE, rm_files = TRUE, asr = TRUE, rep = 1, dir = NULL,
                       n_starts = NULL, ...){
  exec <- path.expand(exec)
  if(file.access(exec, mode=1) == -1) {
    stop("The file ", exec, " cannot be executed.")
  }
  version_test <- grepl("PTHREADS", exec)
  if(version_test){
    stop("Dowser currently only supports the raxmlHPC based options. Please reinitate the function using the 'raxmlHPC' based executable.")
  }
  version_test <- grepl("-ng", exec)
  if(!version_test){
    stop("Please use raxml-ng, not an older version of RAxML.")
  }
  if(!is.null(partition)){
    if(seq != "hlsequence"){
      stop("RAxML partition models are currently only supported for heavy and light chain partitions. Please include clones that contain both.")
    }
  }
  if(!is.null(dir)){
    dir <- path.expand(dir)
    if(!dir.exists(dir)){
      dir.create(dir)
    }
  }else{
    dir <- alakazam::makeTempDir(name)
  }
  clone_seqids <- clone@data$sequence_id
  clone_seqids[length(clone_seqids)+1] <- "Germline"
  if(seq == "hlsequence"){
    clone_seqs <- clone@data$hlsequence
    g <- clone@hlgermline
  }else if(seq == "sequence"){
    clone_seqs <- clone@data$sequence
    g <- clone@germline
  }else if(seq == "lsequence"){
    clone_seqs <- clone@data$lsequence
    g <- clone@lgermline
  }else{
    stop(paste(seq, "not a recognized sequence type"))
  }
  clone_seqs[length(clone_seqs)+1] <- g
  
  # check to make sure that the adjusted sequences are the same length
  nchar_seqs <- nchar(clone_seqs)
  if(nrow(data.frame(table(nchar_seqs))) > 1){
    stop("Sequence lengths do not match.")
  }
  if(length(clone_seqids) != length(clone_seqs)){
    stop("The number of sequences does not match the number of sequence ids.")
  }
  if(!is.numeric(rseed)){
    stop("The random seed needs to be a numeric value")
  }
  
  name <- paste0(name, "_", clone@clone, "_", rep)
  # create the data file for raxml
  fileConn<-file(file.path(dir, paste0(name, "_input_data.phy")))
  writeLines(c(paste0(as.character(length(clone_seqs)), " ", 
                      as.character(nchar(clone_seqs[1])))), fileConn)
  for(i in 1:length(clone_seqs)){
    write(paste0(clone_seqids[i], "    ", clone_seqs[i]), file=file.path(dir, paste0(name, "_input_data.phy")), 
          append = TRUE)
  }
  close(fileConn)
  input_data <- file.path(dir, paste0(name, "_input_data.phy"))
  
  command <- paste("--model", model, "--seed", rseed, "-msa", 
                   input_data, "-prefix", paste0(dir,"/", name), "--threads 1",
                   "--data-type", data_type, "--force msa")
  if(!is.null(n_starts)){
    command <- paste0(command, " --tree pars{",n_starts,"}")
  }
  if(!is.null(starting_tree)){
    if(from_getTrees){
      ape::write.tree(starting_tree, file.path(dir, paste0(name, "_og_starting_tree.tree")))
      starting_tree <- ape::read.tree(file.path(dir, paste0(name, "_og_starting_tree.tree")))
    }
    command <- paste(command, "--tree", starting_tree)
  }
  if(!is.null(partition)){
    heavy_index <- clone@locus == "IGH"
    end_heavy <- paste(strsplit(clone_seqs,split="")[[1]][heavy_index], collapse = "")
    fileConn<-file(file.path(dir, paste0(name, "_partition.txt")))
    write(paste0(model, ", p1 = 1-", nchar(end_heavy)), file=file.path(dir, paste0(name, "_partition.txt")), 
          append = TRUE)
    write(paste0(model, ", p2 = ", nchar(end_heavy)+1, "-", nchar(clone_seqs[1])), 
          file=file.path(dir, paste0(name, "_partition.txt")), 
          append = TRUE)
    close(fileConn)
    old_command <- strsplit(command, "--seed")[[1]][2]
    new_model <- paste("--model", file.path(dir, paste0(name, "_partition.txt")), "--seed")
    command <- paste0(new_model, old_command, " --brlen ", partition)
  }
  
  params <- list(exec, command, stdout=TRUE, stderr=TRUE)
  
  status <- tryCatch(do.call(base::system2, params), error=function(e){
    print(paste("RAxML error, trying again: ",e));
    return(e)
  }, warning=function(w){
    print(paste("RAxML warnings, trying again: ",w));
    return(w)
  })
  
  if(length(status) == 0){
    stop("RAxML was not run. Check input data.")
  }
  
  # check if there is a reduced file -- if so the file creation step failed
  if(file.exists(file.path(dir,paste0(name, "_input_data.phy.reduced")))){
    stop("Preprocessing broke. Not all noninformative sites by RAxML's definition were removed.")
  }
  
  if(asr){
    tree <- rerootTree(ape::unroot(ape::read.tree(file.path(dir,paste0(name, ".raxml.bestTree")))), "Germline", verbose = 0)
    starting_tree <- file.path(dir, paste0(name, "_rerooted.tree"))
    ape::write.tree(tree, starting_tree)
    command <- paste("--model", model, "--seed", rseed, "-msa", 
                     input_data, "-prefix", paste0(dir,"/", name, "_asr"), "--threads 1",
                     "--tree", starting_tree, "--ancestral", "data-type", data_type, 
                     "--force msa")
    if(!is.null(partition)){
      old_command <- strsplit(command, "--seed")[[1]][2]
      new_model <- paste("--model", file.path(dir, paste0(name, "_partition.txt")), "--seed")
      command <- paste0(new_model, old_command, " --brlen ", partition)
    }
    
    params <- list(exec, command, stdout=TRUE, stderr=TRUE)
    status <- tryCatch(do.call(base::system2, params), error=function(e){
      print(paste("RAxML error, trying again: ",e));
      return(e)
    }, warning=function(w){
      print(paste("RAxML warnings, trying again: ",w));
      return(w)
    })
    
    # check that topology is the same
    asr_tree <- ape::read.tree(file.path(dir, paste0(name, "_asr.raxml.ancestralTree")))
    difference_check <- phangorn::RF.dist(asr_tree, tree)
    if(difference_check > 0){
      stop("ASR failed. Retry")
    }
    
    # use the asr_tree from here on out
    tree <- asr_tree
    results <- list()
    results$clone <- clone@clone
    results$nseq <- length(clone@data[[seq]]) 
    results$nsite <- nchar(clone@data[[seq]][1]) 
    results$tree_length <- sum(tree$edge.length)
    if(!is.null(partition) && partition == "unlinked"){
      p_trees <- ape::read.tree(file.path(dir, paste0(name, ".raxml.bestPartitionTrees")))
      p1 <- sum(p_trees[[1]]$edge.length)
      p2 <- sum(p_trees[[2]]$edge.length)
      tree_length <- c(paste("Best Tree:", results$tree_length), paste("Heavy Chain Tree:", p1),
                       paste("Light Chain Tree:", p2))
      results$tree_length <- tree_length
    }
    bestmodel <- readLines(file.path(dir, paste0(name, ".raxml.bestModel")))
    likelihood <- readLines(file.path(dir, paste0(name, "_asr.raxml.log")))
    results$likelihood <- as.numeric(strsplit(likelihood[grep("final logLikelihood:", likelihood)],
                                              "final logLikelihood: ")[[1]][2])
    results$model <- bestmodel
    tree$parameters <- results
    # get the ASR for the nodes
    nnodes <- length(unique(c(tree$edge[,1],tree$edge[,2])))
    asr_seqs <- readLines(file.path(dir, paste0(name, "_asr.raxml.ancestralStates")))
    # find the max node 
    # max_node <- strsplit(asr_seqs, "Node")
    # max_node <- length(max_node)
    # if(nnodes != length(tree$tip.label) + max_node){
    #   stop("Internal node error")
    # }
    tree$nodes <- rep(list(sequence=NULL),times=nnodes)

    #KBH 10/9/24 corrected for ASR matching - please check below
    asr_seqs <- strsplit(asr_seqs, split="\t")
    seqs <- sapply(asr_seqs, function(x)x[2])
    names <- sapply(asr_seqs, function(x)x[1])
    asr_seqs <- seqs
    names(asr_seqs) <- names
    
    # CGJ 10/14/24 add in the internal germline node as "" to node.label and what not
    tree$node.label <- unlist(lapply(0:length(tree$node.label), function(x){
      if(x == 0){
        return("")
      } else{ 
        value <- tree$node.label[x]
      }
    }))
    
    ASR <- list()
    for(i in 1:nnodes){
      if(i <= length(tree$tip.label)){
        # find the ith value in tip.label (assiocated with 1-x on the edge table)
        seq_id <- tree$tip.label[i]
        if(seq_id == "Germline"){
          if(clone@phylo_seq == "sequence"){
            asr_seq <- clone@germline
          } else if(clone@phylo_seq == "hlsequence"){
            asr_seq <- clone@hlgermline
          } else{
            asr_seq <- clone@lgermline
          }
        } else{
          asr_seq <- clone@data[[clone@phylo_seq]][which(clone@data$sequence_id == seq_id)]
        }
        names(asr_seq) <- seq_id
      } else{
        label <- tree$node.label[i - length(tree$tip.label)]
        if(label == ""){
          if(clone@phylo_seq == "sequence"){
            asr_seq <- clone@germline
          } else if(clone@phylo_seq == "hlsequence"){
            asr_seq <- clone@hlgermline
          } else{
            asr_seq <- clone@lgermline
          }
        } else{
          asr_seq <- asr_seqs[label]
          asr_seq <- gsub("-", "N", asr_seq)
        }
      }
      ASR[[i]] <- asr_seq
    }

    
    tree$nodes <- lapply(1:length(tree$nodes),function(x){
      tree$nodes[[x]]$sequence <- ASR[[x]]
      tree$nodes[[x]]
    })
    
    tree <- rerootTree(tree, "Germline", verbose=0)
  }else {
    tree <- rerootTree(ape::unroot(ape::read.tree(file.path(dir,paste0(name, ".raxml.bestTree")))), "Germline", verbose = 0)
    tree$germid <- paste0(clone@clone, "_GERM")
    results <- list()
    results$clone <- clone@clone
    results$nseq <- length(clone@data[[seq]]) 
    results$nsite <- nchar(clone@data[[seq]][1])
    results$tree_length <- sum(tree$edge.length) 
    nnodes <- length(unique(c(tree$edge[,1],tree$edge[,2])))
    tree$nodes <- rep(list(sequence=NULL),times=nnodes)
    if(!is.null(partition) && partition == "unlinked"){
      p_trees <- ape::read.tree(file.path(dir, paste0(name, ".raxml.bestPartitionTrees")))
      p1 <- sum(p_trees[[1]]$edge.length)
      p2 <- sum(p_trees[[2]]$edge.length)
      tree_length <- c(paste("Best Tree:", results$tree_length), paste("Heavy Chain Tree:", p1),
                       paste("Light Chain Tree:", p2))
      results$tree_length <- tree_length
      for(i in 1:length(tree$tip.label)){
        tip <- tree$tip.label[i]
        if(tip != "Germline"){
          seq_num <- which(clone@data$sequence_id == tip)
          if(seq == "hlsequence"){
            tree$nodes[i]$sequence <- clone@data$hlsequence[seq_num]
          } else if(seq == "sequence"){
            tree$nodes[i]$sequence <- clone@data$sequence[seq_num]
          }else{
            tree$nodes[i]$sequence <- clone@data$lsequence[seq_num]
          }
        }else{
          if(seq == "hlsequence"){
            tree$nodes[i]$sequence <- clone@hlgermline
          } else if(seq == "sequence"){
            tree$nodes[i]$sequence <- clone@germline
          }else{
            tree$nodes[i]$sequence <- clone@lgermline
          }
        }
      }
    }
    bestmodel <- readLines(file.path(dir, paste0(name, ".raxml.bestModel")))
    likelihood <- readLines(file.path(dir, paste0(name, ".raxml.log")))
    results$likelihood <- as.numeric(strsplit(likelihood[grep("Final LogLikelihood:", likelihood)],
                                              "Final LogLikelihood: ")[[1]][2])
    results$model <- bestmodel
    tree$parameters <- results
    tipseqs <- clone@data[[seq]]
    names(tipseqs) <- clone@data$sequence_id
    if(seq == "sequence"){
      tipseqs <- c(tipseqs,"Germline"=clone@germline)
    }else if(seq == "lsequence"){
      tipseqs <- c(tipseqs,"Germline"=clone@lgermline)
    }else if(seq == "hlsequence"){
      tipseqs <- c(tipseqs,"Germline"=clone@hlgermline)
    }else{
      stop(paste("phylo_seq not recognized",clone))
    }
    if(sum(!tree$tip.label %in% names(tipseqs)) != 0 ||
       sum(!names(tipseqs) %in% tree$tip.label) != 0){
      stop(paste("Tip sequences do not match in clone",clone@clone))
    }
    for(i in 1:nnodes){
      if(i <= length(tree$tip.label)){
        tree$nodes[i]$sequence <- tipseqs[tree$tip.label[i]]
      }
    }
  }
  tree$name <- clone@clone
  tree$seq <- clone@phylo_seq
  tree_method <- readLines(file.path(dir, paste0(name, ".raxml.log")))
  tree$tree_method <- paste0("RAxML:", strsplit(tree_method[grep("Model: ", tree_method)],
                               "Model: ")[[1]][2])
  tree$edge_type <- "genetic_distance"
  if(rm_files){
    unlink(file.path(dir, paste0(name,"*")))
  }
  return(tree)
}


#' Reroot phylogenetic tree to have its germline sequence at a zero-length branch
#' to a node which is the direct ancestor of the tree's UCA. Assigns \code{uca}
#' to be the ancestral node to the tree's germline sequence, as \code{germid} as
#' the tree's germline sequence ID.
#'
#' @param   tree      An ape \code{phylo} object
#' @param   germline  ID of the tree's predicted germline sequence
#' @param   min       Maximum allowed branch length from germline to root
#' @param   verbose    amount of rubbish to print
#' @return  \code{phylo} object rooted at the specified germline
#' 
#' @export
rerootTree <- function(tree, germline, min=0.001, verbose=1){
  ntip <- length(tree$tip.label)
  uca <- ntip+1
  if(!germline %in% tree$tip.label){
    stop(paste(germline,"not found in tip labels!"))
  }
  olength <- sum(tree$edge.length)
  odiv <- ape::cophenetic.phylo(tree)[germline,]
  if(ape::is.rooted(tree)){
    if(verbose > 0){
      print("unrooting tree!")
    }
    root <- ape::getMRCA(tree,
                         tip=tree$tip.label)
    max <- max(tree$edge)
    rindex <- tree$edge[,1] == root
    redge <- tree$edge[rindex,]
    parent <- redge[1,2]
    child <- redge[2,2]
    if(parent <= length(tree$tip.label)){
      parent <- redge[2,2]
      child <- redge[1,2]
    }
    if(tree$tip.label[child] == germline &&
       tree$edge.length[tree$edge[,1] == root &
                        tree$edge[,2] == child] <= min){
      if(verbose > 0){
        print("tree already rooted at germline!")
      }
      return(tree)
    }
    warning("Rooting already rooted trees not fully supported.")
    # if tree rooted somewhere besides the germline
    # cut out the root node and directly attach
    # its former descendants to each other.
    tree$edge <- tree$edge[!rindex,]
    tree$edge <- rbind(tree$edge,c(parent,child))
    sumedge <- sum(tree$edge.length[rindex])
    tree$edge.length <- tree$edge.length[!rindex]
    tree$edge.length[length(tree$edge.length)+1] <- sumedge
    tree$Nnode <- tree$Nnode - 1
    if(parent != uca){
      #if new parent doesn't have correct uca number, swap it
      #uca node should have been the node that was cut out
      #if not the parent
      if(uca %in% tree$edge){
        stop("something weird happened during unrooting")
      }
      root <- parent
      tree$edge[tree$edge[,1]==parent,1] <- uca
      tree$edge[tree$edge[,2]==parent,2] <- uca
      if(!is.null(tree$nodes)){
        s <- tree$nodes[[uca]]
        tree$nodes[[uca]] <- tree$nodes[[parent]]
        tree$nodes[[parent]] <- s
      }
    }
    # set the node information for the appropriate root
    tree$edge[tree$edge[,1]==max,1] <- root
    tree$edge[tree$edge[,2]==max,2] <- root
    if(!is.null(tree$nodes)){
      tree$nodes[[max]] <- tree$nodes[[root]]
      tree$nodes[[max]] <- NULL
    }
  }
  #pairwise patristic distance among all nodes
  odist <- ape::dist.nodes(tree)
  
  edge <- tree$edge
  germid <- which(tree$tip.label == germline)
  
  max <- max(edge)
  nnode <- max+1 #new node number
  uca <- ntip+1 #uca needs to be first internal node
  
  # replace current uca (hereafter mrca) with new node number in edge list
  edge[edge[,1]==uca,1] <- nnode
  edge[edge[,2]==uca,2] <- nnode
  
  # make MRCA connect to new UCA node instead of germline
  edge[edge[,2] == germid,2] <- uca
  
  #add 0 length edge from UCA to germline
  edge <- rbind(edge,c(uca,germid)) 
  tree$edge.length[length(tree$edge.length)+1] <- 0

  # recursive function to swap nodes while
  # preserving internal node ids
  # this is what actually reroots the trees
  swap <- function(tnode, edge, checked){
    if(tnode %in% checked){
      print("r edge")
      return(edge)
    }
    checked <- c(checked,tnode)
    children <- edge[edge[,1] == tnode,2] #get children of this node
    parent <- edge[edge[,2] == tnode,1] #get parent of this node
    # only one child node, or if parent hasn't been checked, swap target and parent
    if(length(children) < 2 || sum(!parent %in% checked) > 0){
      parent <- edge[edge[,2] == tnode,1]
      parent <- parent[!parent %in% checked]
      edge[edge[,1] == parent & edge[,2] == tnode,] <- c(tnode,parent)
      children <- edge[edge[,1] == tnode,2]
    }
    #make sure descendant branches are facing the correct direction
    for(tnode in children){
      if(!tnode %in% checked){
        edge <- swap(tnode,edge,checked)
      }
    }
    return(edge)
  }
  
  #place new UCA node at root of tree
  edge <- swap(uca, edge, checked=c(1:ntip))
  
  tree$edge <- edge
  tree$Nnode <- length(unique(edge[,1]))
  tree <- ape::reorder.phylo(tree,"postorder")
  
  # swap node information between mrca and old uca
  # uca gets same info as germline
  if(!is.null(tree$nodes)){
    tree$nodes[[nnode]] <- tree$nodes[[uca]]
    tree$nodes[[uca]] <- tree$nodes[[germid]]
  }
  
  # reset and re-order tree
  attr(tree, "order") <- NULL
  tree <- ape::ladderize(tree, right=FALSE)
  
  # sanity check tree length, divergence, and internal node distances
  nlength <- sum(tree$edge.length)
  ndiv <- ape::cophenetic.phylo(tree)[germline,]
  if(abs(nlength - olength) > 0.001){
    stop(paste("Error in rerooting tree",tree$name,
               "tree length not consistent"))
  }
  divergence_diff <- odiv - ndiv[names(odiv)]
  if(max(abs(divergence_diff)) > 0.001){
    stop(paste("Error in rerooting tree",tree$name,
               "germline divergences not consistent"))
  }
  ndist <- ape::dist.nodes(tree)
  ndist[uca,] <- ndist[nnode,]
  ndist[,uca] <- ndist[,nnode]
  ndist <- ndist[,-nnode]
  ndist <- ndist[-nnode,]
  max_diff <- max(odist - ndist)
  if(abs(max_diff) > 0.001){
    stop(paste("Error in rerooting tree",tree$name,
               "node distances not consistent"))
  }
  return(tree)
}

#' Estimate lineage tree topologies, branch lengths,
#' and internal node states if desired
#'
#' \code{getTrees} Tree building function.
#' @param    clones     a tibble of \code{airrClone} objects, the output of
#'                      \link{formatClones}
#' @param    trait      trait to use for parsimony models (required if
#'                      \code{igphyml} specified)
#' @param    build      program to use for tree building (pratchet, pml,
#'                      dnapars, dnaml, igphyml, raxml)
#' @param    exec       location of desired phylogenetic executable
#' @param    igphyml    optional location of igphyml executible for parsimony
#' @param    id         unique identifer for this analysis (required if
#'                      \code{igphyml} or \code{dnapars} specified)
#' @param    dir        directory where temporary files will be placed.
#' @param    modelfile  file specifying parsimony model to use
#' @param    fixtrees   if TRUE, use supplied tree topologies
#' @param    nproc      number of cores to parallelize computations
#' @param    quiet      amount of rubbish to print to console
#' @param    rm_temp    remove temporary files (default=TRUE)
#' @param    palette    deprecated
#' @param    seq        column name containing sequence information
#' @param    collapse   Collapse internal nodes with identical sequences?
#' @param    ...        Additional arguments passed to tree building programs
#'
#' @return   A list of \code{phylo} objects in the same order as \code{data}.
#'
#' @details
#' Estimates phylogenetic tree topologies and branch lengths for a list of
#' \code{airrClone} objects. By default, it will use phangnorn::pratchet to
#' estimate maximum parsimony tree topologies, and ape::acctran to estimate
#' branch lengths. If \code{igpyhml} is specified, internal node \code{trait}
#' values will be predicted by maximum parsimony. In this case, \code{dir} will
#' need to be specified as a temporary directory to place all the intermediate
#' files (will be created if not available). Further, \code{id} will need to
#' specified to serve as a unique identifier for the temporary files. This
#' should be chosen to ensure that multiple \code{getTrees} calls using the same
#' \code{dir} do not overwrite each others files.
#'
#' \code{modelfile} is written automatically if not specified, but doesn't
#' include any constraints. Intermediate files are deleted by default. This can
#' be toggled using (\code{rm_files}).
#'
#' For examples and vignettes, see https://dowser.readthedocs.io
#'
#' @seealso \link{formatClones}, \link{findSwitches}, \link{buildPhylo},
#' \link{buildPratchet}, \link{buildPML}, \link{buildIgphyml}, \link{buildRAxML}
#' @examples
#' data(ExampleClones)
#' trees <- getTrees(ExampleClones[10,])
#' plotTrees(trees)[[1]]
#'
#' \dontrun{
#' data(ExampleClones)
#'
#' trees <- getTrees(ExampleClones[10,],igphyml="/path/to/igphyml",
#'          id="temp",dir="temp", trait="sample_id")
#' plotTrees(trees)[[1]]
#' }
#' @export
getTrees <- function(clones, trait=NULL, id=NULL, dir=NULL,
                     modelfile=NULL, build="pratchet", exec=NULL, igphyml=NULL,
                     fixtrees=FALSE, nproc=1, quiet=0, rm_temp=TRUE, palette=NULL,
                     seq=NULL, collapse=FALSE, ...){

  if(is.null(exec) && (!build %in% c("pratchet", "pml"))){
    stop("exec must be specified for this build option")
  }
  if(!is.null(dir)){
    dir <- path.expand(dir)
  }
  data <- clones$data
  if(!inherits(data, "list")){
    data <- list(data)
  }
  if(!inherits(data[[1]], "airrClone")){
    stop("Input data must be a list of airrClone objects")
  }
  if(!is.null(palette)){
    warning("palette option is deprecated in getTrees, specify in plotTrees")
    palette <- NULL
  }
  
  # make sure all sequences and germlines within a clone are the same length
  unlist(lapply(data, function(x){
    if(x@phylo_seq == "hlsequence"){
      germline <- x@hlgermline
      seqs <- x@data$hlsequence
    }else if(x@phylo_seq == "lsequence"){
      germline <- x@lgermline
      seqs <- x@data$lsequence
    }else{
      germline <- x@germline
      seqs <- x@data$sequence
    }
    if(any(nchar(germline) != nchar(seqs))){
      stop(paste0("Sequence and/or germline lengths of clone ",
                  x@clone," are not equal."))
    }
  }))

  if(fixtrees){
    if(!"trees" %in% names(clones)){
      stop("trees column must be specified if fixtrees=TRUE")
    }
    if(!inherits(clones$trees[[1]], "phylo")){
      stop("Trees must be a list of class phylo")
    }
    trees <- clones$trees
  }else{
    trees <- NULL
  }

  if(is.null(id)){
    id <- "sample"
  }
  big <- FALSE
  if(sum(unlist(lapply(data, function(x)nrow(x@data)))) > 10000){
    big <- TRUE
  }
  if(!rm_temp && big){
    warning("Large dataset - best to set rm_temp=TRUE")
  }
  if(!is.null(igphyml)){
    igphyml <- path.expand(igphyml)
    if(!is.null(dir)){
      if(!dir.exists(dir)){
        dir.create(dir)
      }
    }else{
      dir <- alakazam::makeTempDir(id)
      if(big){
        warning("Large dataset - best to set dir and id params")
      }
    }
    if(is.null(trait)){
      stop("trait must be specified when igphyml-based trait reconstruction")
    }
    # remove problematic characters from trait values
    data <- lapply(data,function(x){
      x@data[[trait]] <- gsub(":|;|,|=| |_","-",x@data[[trait]])
      x})
    if(is.null(modelfile)){
      states <- unique(unlist(lapply(data,function(x)x@data[,trait])))
      modelfile <- makeModelFile(states,
                                 file=file.path(dir,paste0(id,"_modelfile.txt")))
    }else{
      states <- readModelFile(modelfile)
    }
    #if igphyml is specified, append trait value to sequence ids
    if(!is.null(trees)){
      indexes <- 1:length(data)
      trees <- lapply(indexes,function(x){
        tree <- trees[[x]]
        datat <- data[[x]]
        for(id in datat@data$sequence_id){
          trait_temp <- filter(datat@data,
                               !!rlang::sym("sequence_id")==id)[[trait]]
          tree$tip.label[tree$tip.label == id] <-
            paste0(id,"_",trait_temp)
        }
        tree})
    }
    data <- lapply(data,function(x){
      x@data$sequence_id <- paste0(x@data$sequence_id,"_",x@data[[trait]])
      x})
  }

  if(build=="dnapars" || build=="igphyml" || build=="dnaml" || !is.null(igphyml) || build=="raxml"){
    if(!is.null(dir)){
      if(!dir.exists(dir)){
        dir.create(dir)
      }
    }else{
      dir <- alakazam::makeTempDir(id)
      if(big){
        warning("Large dataset - best to set dir and id params")
      }
    }
  }

  if(!inherits(data, "list")){
    data <- list(data)
  }
  if(!is.null(dir)){
    if(!dir.exists(dir)){
      dir.create(dir)
    }
  }
  rm_dir <- NULL
  if(rm_temp){
    rm_dir=file.path(dir,paste0(id,"_recon_trees"))
  }

  reps <- as.list(1:length(data))

  if(is.null(seq)){
    seqs <- unlist(lapply(data,function(x)x@phylo_seq))
  }else{
    seqs <- rep(seq,length=length(data))
  }
  
  if(build=="dnapars" || build=="dnaml"){
    trees <- parallel::mclapply(reps,function(x)
      tryCatch(buildPhylo(data[[x]],
                          exec=exec,
                          temp_path=file.path(dir,paste0(id,"_trees_",x)),
                          rm_temp=rm_temp,seq=seqs[x],tree=trees[[x]],...),
               error=function(e)e),
      mc.cores=nproc)
  }else if(build=="pratchet"){
    trees <- parallel::mclapply(reps,function(x)
      tryCatch(buildPratchet(data[[x]],seq=seqs[x],
                             tree=trees[[x]],...),error=function(e)e),
      mc.cores=nproc)
  }else if(build=="pml"){
    trees <- parallel::mclapply(reps,function(x)
      tryCatch(buildPML(data[[x]],seq=seqs[x],
                        tree=trees[[x]],...),error=function(e)e),
      mc.cores=nproc)
  } else if(build=="igphyml"){
    if(rm_temp){
      rm_dir <- file.path(dir,id)
    }else{
      rm_dir <- NULL
    }
    trees <-
      tryCatch(buildIgphyml(data,
                            igphyml=exec,
                            temp_path=file.path(dir,id),
                            rm_files=rm_temp,
                            rm_dir=rm_dir,
                            trees=trees,nproc=nproc,id=id,...),error=function(e)e)
    if(inherits(trees, "error")){
      stop(trees)
    }
  } else if(build=="raxml"){ # CGJ 2/20/23
    trees <- parallel::mclapply(reps,function(x)
      tryCatch(buildRAxML(data[[x]], 
                          seq=seqs[x],
                          exec = exec,
                          rm_files = rm_temp,
                          dir=dir,
                          starting_tree=trees[[x]],...),error=function(e)e),
      mc.cores=nproc)
  } else{
    stop("build specification ", build, " not recognized")
  }

  # save points for data have been saved as the following
  # trees -> trees_save
  # data -> data_save
  # clones -> clones_save
  #catch any tree inference errors
  if(build != "igphyml"){
    errors <- unlist(lapply(trees, function(x) inherits(x, "error")))
    messages <- trees[errors]
    errorclones <- clones$clone_id[errors]
    trees <- trees[!errors]
    data <- data[!errors]
    clones <- clones[!errors,]
    if(length(errorclones) > 0){
      warning(paste("Tree building failed for clones",
                    paste(errorclones,collapse=", ")))
      me <- lapply(messages, function(x)warning(x$message))
    }
  }
  if(length(trees) == 0){
    stop("No trees left!")
  }

  # make sure trees, data, and clone objects are in same order
  tree_names <- unlist(lapply(trees, function(x)x$name))
  data_names <- unlist(lapply(data, function(x)x@clone))
  m <- match(tree_names, data_names)
  data <- data[m]
  m <- match(tree_names, clones$clone_id)
  clones <- clones[m,]

  if(build == "igphyml" | build == "raxml"){
    clones$parameters <- lapply(trees,function(x)x$parameters)
  }

  if(!is.null(igphyml)){
    file <- writeLineageFile(data=data, trees=trees, dir=dir,
                             id=id, trait=trait, rep="trees")

    mtrees <- reconIgPhyML(file, modelfile, igphyml=igphyml, 
                           mode="trees", id=NULL, quiet=quiet, nproc=nproc,
                           rm_files=rm_temp, rm_dir=rm_dir, states=states, 
                           palette=palette, ...)
    
    # remove trait value from tips
    mtrees <- lapply(mtrees,function(x){
      ids <- strsplit(x$tip.label,split="_")
      x$tip.label <- unlist(lapply(ids,function(t)
        paste(t[1:(length(t)-1)],collapse="_")))
      x$tip.label[x$tip.label == x$name] <- "Germline"
      x
    })
    
  }else{
    mtrees <- trees
  }
  # Sanity checks
  match <- unlist(lapply(1:length(data), function(x){
    data[[x]]@clone == mtrees[[x]]$name
  }))
  if(sum(!match) > 0){
    stop("Clone and tree names not in proper order!")
  }
  clones$trees <- mtrees
  if(collapse){
    clones <- collapseNodes(clones)
  }
  if(sum(clones$clone_id != 
         unlist(lapply(trees,function(x)x$name))) > 0){
    stop("Tree column names don't match clone IDs")
  }
  
  clones
}



#' Scale branch lengths to represent either mutations or mutations per site.
#' 
#' \code{scaleBranches} Branch length scaling function.
#' @param    clones      a tibble of \code{airrClone} and \code{phylo} objects,
#'                      the output of \link{getTrees}.
#' @param    edge_type  Either \code{genetic_distance} (mutations per site) or 
#'                      \code{mutations}
#' 
#' @return   A tibble with \code{phylo} objects that have had branch lengths 
#'           rescaled as specified.
#'
#' @details
#' Uses clones$trees[[1]]$edge_type to determine how branches are currently scaled.
#'  
#' @seealso \link{getTrees}
#' @export
scaleBranches <- function(clones, edge_type="mutations"){
  if(!"tbl" %in% class(clones)){
    print(paste("clones is of class",class(clones)))
    stop("clones must be a tibble of airrClone objects!")
  }else{
    if(!inherits(clones$data[[1]], "airrClone")){
      print(paste("clones is list of class",class(clones$data[[1]])))
      stop("clones must be a list of airrClone objects!")
    }
  }
  if(!"trees" %in% names(clones)){
    stop("clones must have trees column!")
  }
  lengths <- unlist(lapply(1:length(clones$trees),
                           function(x){
                             if(clones$data[[x]]@phylo_seq == "hlsequence"){
                               return(nchar(clones$data[[x]]@hlgermline))
                             }else if(clones$data[[x]]@phylo_seq == "lsequence"){
                               return(nchar(clones$data[[x]]@lgermline))
                             }else{
                               return(nchar(clones$data[[x]]@germline))
                             }}))
  
  trees <- lapply(1:length(clones$trees),function(x){
    if(clones$trees[[x]]$edge_type == "mutations" && 
       edge_type == "genetic_distance"){
      clones$trees[[x]]$edge.length <- 
        clones$trees[[x]]$edge.length/lengths[x]
      clones$trees[[x]]$edge_type <- "genetic_distance"
      clones$trees[[x]]
    }else if(clones$trees[[x]]$edge_type == "genetic_distance" &&
             edge_type == "mutations"){
      clones$trees[[x]]$edge.length <- 
        clones$trees[[x]]$edge.length*lengths[x]
      clones$trees[[x]]$edge_type <- "mutations"
      clones$trees[[x]]
        }else if(clones$trees[[x]]$edge_type == "genetic_distance_codon" &&
             edge_type == "mutations"){
      clones$trees[[x]]$edge.length <-
        clones$trees[[x]]$edge.length*lengths[x]/3
      clones$trees[[x]]$edge_type <- "mutations"
      clones$trees[[x]]
    }else{
      warning("edge conversion type not yet supported")
      clones$trees[[x]]
    }})
  
  clones$trees <- trees
  clones
}

#' Collapse internal nodes with the same predicted sequence
#' 
#' \code{collapseNodes} Node collapsing function.
#' @param    trees    a tibble of \code{airrClone} objects, the output of \link{getTrees}
#' @param    tips     collapse tips to internal nodes? (experimental)  
#' @param    check    check that collapsed nodes are consistent with original tree
#' 
#' @return   A tibble with \code{phylo} objects that have had internal nodes collapsed.
#'
#' @details
#' Use plotTrees(trees)[[1]] + geom_label(aes(label=node)) + geom_tippoint() to show
#' node labels, and getSeq to return internal node sequences
#'  
#' @seealso \link{getTrees}
#' @export
collapseNodes <- function(trees, tips=FALSE, check=TRUE){
  if(!"phylo" %in% class(trees)){
    trees$trees <- lapply(trees$trees,function(x)
      collapseNodes(x,tips))
    return(trees)
  }
  otrees <- trees
  edges <- trees$edge
  ttips <- trees$tip.label
  edge_l <- trees$edge.length
  btip <- length(trees$tip.label)
  if(is.null(trees$node.label)){
    trees$node.label <- rep("",length(unique(edges)))
  }
  maxiter <- 1000
  while(maxiter > 0){
    # edge list, plus whether or not reconstructed sequences are identical
    edges <- cbind(edges[,1:2],unlist(lapply(1:nrow(edges),function(x)
      trees$nodes[[edges[x,1]]]$sequence==trees$nodes[[edges[x,2]]]$sequence)))
    # get list of edges to collapse
    if(!tips){
      zedges <- edges[edges[,2] > btip & edges[,3] == 1,]
    }else{
      mrca <- ape::getMRCA(trees,tip=trees$tip.label)
      zedges <- edges[edges[,1] != mrca & edges[,3] == 1,]
    }
    if(length(zedges) > 0){
      if(!is.matrix(zedges)){
        dim(zedges) <- c(1,3)
      }
    }else{
      break
    }
    # replace identical nodes with smaller of the node numbers
    ms <- c() #smaller node numbers
    ns <- c() #larger node numbers
    for(i in 1:nrow(zedges)){
      m <- min(zedges[i,1:2])
      n <- max(zedges[i,1:2])
      edges[edges[,1] == n,1] <- m
      edges[edges[,2] == n,2] <- m
      ms <- c(ms,m)
      ns <- c(ns,n)
    }
    # remove edges that lead to their own node
    edge_l <- edge_l[edges[,1] != edges[,2]]
    edges <- edges[edges[,1] != edges[,2],]
    maxiter <- maxiter - 1
  }
  if(maxiter == 0){
    stop("Exceeded maximum node collapse iterations")
  }
  # make collapsed tips into internal nodes
  if(tips){
    tnodes <- unique(edges[edges[,1] <= btip,1])
    m <- max(edges)+1
    for(n in tnodes){ #replace tip number with new number
      edges[edges[,1] == n,1] <- m
      edges[edges[,2] == n,2] <- m
      trees$nodes[[m]] <- trees$nodes[[n]]
      trees$nodes[[m]]$id <- ttips[n]
      m <- m + 1
    }
    ttips <- ttips[-tnodes]
  }
  # make key of new node names
  nodes <- sort(unique(c(edges[,1],edges[,2])))
  nnodes <- 1:length(nodes)
  names(nnodes) <- nodes
  
  # replace old nodes with new node names and info
  nedges <- edges
  nsequences <- as.list(1:length(nodes))
  for(n in names(nnodes)){
    new <- nnodes[n]
    nedges[edges[,1] == as.numeric(n),1] <- new
    nedges[edges[,2] == as.numeric(n),2] <- new
    nsequences[[new]] <- trees$nodes[[as.numeric(n)]]
    if(is.null(nsequences[[new]]$id)){
      nsequences[[new]]$id <- ""
    }
  }
  trees$edge <- nedges[,1:2]
  trees$edge.length <- edge_l
  trees$nodes <- nsequences
  trees$Nnode <- length(nsequences)-length(ttips)
  trees$tip.label <- ttips
  nodelabs <- unlist(lapply(trees$nodes,function(x){
    if(is.null(x$id)){
      ""
    }else{
      x$id
    }}))
  trees$node.label <- nodelabs[(length(ttips)+1):length(nsequences)]
  
  # reset and re-order tree
  attr(trees, "order") <- NULL
  trees <- ape::ladderize(trees, right=FALSE)
  
  if(check){
    if(tips){
      warning("Cannot check nodes while collapsing tips (tips=TRUE)")
    }else{
      subt <- ape::subtrees(otrees)
      for(sub in subt){
        seqs <- sub$tip.label
        node <-  ape::getMRCA(otrees,tip=seqs)
        cnode <- ape::getMRCA(trees,tip=seqs)
        oseq <- otrees$nodes[[node]]$sequence
        nseq <- trees$nodes[[cnode]]$sequence
        if(oseq != nseq){
          stop(paste("Node",node,"in clone tree",trees$name,
                     "does not equal corresponding sequence in collapsed tree"))
        }
      }
    }
  }
  return(trees)
}

#' Return IMGT gapped sequence of specified tree node
#' 
#' \code{getNodeSeq} Sequence retrieval function.
#' @param    data    a tibble of \code{airrClone} objects, the output of 
#'                   \link{getTrees}
#' @param    node    numeric node in tree (see details)
#' @param    tree    a \code{phylo} tree object containing \code{node}
#' @param    clone   if \code{tree} not specified, supply clone ID in \code{data}
#' @param    gaps    add IMGT gaps to output sequences?
#' @return   A vector with sequence for each locus at a specified \code{node}
#'           in \code{tree}.
#'
#' @details
#' Use plotTrees(trees)[[1]] + geom_label(aes(label=node))+geom_tippoint() to show
#' node labels, and getNodeSeq to return internal node sequences
#'  
#' @seealso \link{getTrees}
#' @export
getNodeSeq <- function(data, node, tree=NULL, clone=NULL, gaps=TRUE){
  if(is.null(tree)){
    if(is.null(clone)){
      stop("must provide either tree object or clone ID")
    }
    tree <- dplyr::filter(data,!!rlang::sym("clone_id")==clone)$trees[[1]]
  }
  clone <- dplyr::filter(data,!!rlang::sym("clone_id")==tree$name)$data[[1]]
  seqs <- c()
  seq <- strsplit(tree$nodes[[node]]$sequence,split="")[[1]]
  loci <- unique(clone@locus)
  for(locus in loci){
    if(length(seq) < length(clone@locus)){
      warning("Sequences are shorter than chain vector. Exiting")
    }
    if(length(seq) > length(clone@locus)){
      stop("Sequences are longer than chain vector. Exiting")
    }
    lseq <- seq[clone@locus == locus]
    lseq[is.na(lseq)] <- "N"
    if(gaps){
      nums <- clone@numbers[clone@locus == locus]
      nseq <- rep(".",max(nums))
      nseq[nums] <- lseq
      lseq <- nseq
    }
    seqs <- c(seqs,paste(lseq,collapse=""))
  }
  names(seqs) <- loci
  return(seqs)
}

#' Deprecated! Use getNodeSeq
#' 
#' \code{getSeq} Sequence retrieval function.
#' @param    data    a tibble of \code{airrClone} objects, the output of
#'                   \link{getTrees}
#' @param    node    numeric node in tree (see details)
#' @param    tree    a \code{phylo} tree object containing \code{node}
#' @param    clone   if \code{tree} not specified, supply clone ID in \code{data}
#' @param    gaps    add IMGT gaps to output sequences?
#' @return   A vector with sequence for each locus at a specified \code{node}
#'           in \code{tree}.
#'
#' @seealso \link{getTrees}
#' @export
getSeq <- function(data, node, tree=NULL, clone=NULL, gaps=TRUE){
  
  warning("getSeq is depracated and will be removed. Use getNodeSeq instead.")
  
  return(getNodeSeq(data=data, node=node, tree=tree, clone=clone, gaps=gaps))
}

#' \code{downsampleClone} Down-sample clone to maximum tip/switch ratio
# TODO: Add support for weighting by collapseCount?
#' @param    clone       an \link{airrClone} object
#' @param    trait       trait considered for rarefaction
#'                       \link{getTrees}
#' @param    tree        a \code{phylo} tree object correspond to \code{clone}
#' @param    tip_switch  maximum tip/switch ratio
#' @return   A vector with sequence for each locus at a specified \code{node}
#'           in \code{tree}.
#' @export
downsampleClone <- function(clone, trait, tip_switch=20, tree=NULL){
  
  cdata <- clone@data
  if(!trait %in% names(cdata)){
    stop(paste(trait,"not found in clone data columns"))
  }
  states <- unique(cdata[[trait]])
  
  if(sum(is.na(states) > 0)){
    stop("NA trait values detected, must be removed before trait analysis.")
  }
  
  if(length(states) > 1){
    # if at least one of each state is preserved, there is a minimum
    # of length(states)-1 switches, and a minimum of length(states) tips
    if(tip_switch < length(states)/(length(states)-1)){
      stop(paste("clone",clone@clone,
                 "tip/switch ratio =",tip_switch,"not possible with",
                 length(states),"states"))
    }
    
    ntips <- nrow(cdata)
    
    # randomly select one sequence of each type
    saved <- unlist(lapply(states, function(x)
      sample(cdata[cdata[[trait]] == x,]$sequence_id,size=1)))
    
    # if too many tips, downsample
    if(ntips/(length(states)-1) > tip_switch){
      target <- tip_switch*(length(states)-1)
      tips <- cdata$sequence_id
      tips <- tips[!tips %in% saved]
      rm <- sample(tips, size=ntips-target)
      cdata <- cdata[!cdata$sequence_id %in% rm,]
      
      # check that results are good
      if(nrow(cdata)/(dplyr::n_distinct(cdata[[trait]])-1) != tip_switch){
        stop(paste("clone",clone@clone,"downsampling failed!"))
      }
      # if tree provided, drop selected tips
      if(!is.null(tree)){
        od <- getDivergence(tree)
        # CGJ 9/5/23 found via checkUsagePackage("dowser")
        # otree <- tree
        tree <- ape::drop.tip(tree, tip=rm)
        nd <- getDivergence(tree)
        maxdiff <- max(nd - od[names(nd)])
        meandiff <- mean(nd - od[names(nd)])
        if(maxdiff > 0.001){
          badtip <- names(which.max(nd - od[names(nd)]))
          warning(paste("clone",clone@clone,
                        "downsampling divergences differ by max",
                        signif(maxdiff,digits=2),"mean",
                        signif(meandiff,digits=2), "worst tip", badtip))
        }
        if(sum(!cdata$sequence_id %in% tree$tip.label) == 0 &
           sum(!tree$tip.label %in% c(cdata$sequence_id, "Germline"))){
          stop(paste("clone",clone@clone," tree downsampling failed!")) 
        }
      }
    }
  }
  clone@data <- cdata
  results <- list()
  results$clone <- clone
  results$tree <- tree
  results
}

#' Create a bootstrap distribution for clone sequence alignments, and estimate 
#' trees for each bootstrap replicate.
#' 
#' \code{findSwitches} Phylogenetic bootstrap function.
#' @param clones         tibble \code{airrClone} objects, the output of 
#'                      \link{formatClones}
#' @param permutations    number of bootstrap replicates to perform
#' @param trait         trait to use for parsimony models
#' @param igphyml       location of igphyml executible 
#' @param build         program to use for tree building (phangorn, dnapars)
#' @param exec          location of desired phylogenetic executable
#' @param id            unique identifer for this analysis (required if 
#'                      \code{igphyml} or \code{dnapars} specified)
#' @param dir           directory where temporary files will be placed (required
#'                      if \code{igphyml} or \code{dnapars} specified)
#' @param modelfile     file specifying parsimony model to use
#' @param fixtrees      keep tree topologies fixed?
#'                      (bootstrapping will not be perfomed)
#' @param nproc            number of cores to parallelize computations
#' @param quiet           amount of rubbish to print to console
#' @param rm_temp       remove temporary files (default=TRUE)
#' @param palette       deprecated
#' @param resolve       how should polytomies be resolved? 
#'                       0=none, 1=max parsminy, 2=max ambiguity + polytomy skipping,
#'                       3=max ambiguity
#' @param keeptrees     keep trees estimated from bootstrap replicates? (TRUE)
#' @param lfile         lineage file input to igphyml if desired (experimental)
#' @param rep           current bootstrap replicate (experimental)
#' @param seq           column name containing sequence information
#' @param downsample    downsample clones to have a maximum specified tip/switch ratio?
#' @param tip_switch    maximum allowed tip/switch ratio if downsample=TRUE
#' @param boot_part     is  "locus" bootstrap columns for each locus separately
#' @param force_resolve continue even if polytomy resolution fails?
#' @param ...        additional arguments to be passed to tree building program
#'
#' @return   A list of trees and/or switch counts for each bootstrap replicate.
#'
#' @details
#' Tree building details are the same as \link{getTrees}. 
#' If \code{keeptrees=TRUE} (default) the returned object will contain a list 
#' named "trees" which contains a list of estimated tree objects for each 
#' bootstrap replicate. The object is structured like: 
#' trees[[<replicate>]][[<tree index>]]. If \code{igphyml} is specified 
#' (as well as \code{trait}), the returned object 
#' will contain a \code{tibble} named "switches" containing switch count 
#' information. This object can be passed to \link{testSP} and other functions 
#' to perform parsimony based trait value tests. 
#'
#' Trait values cannot contain values N, UCA, or NTIP. These are reserved for
#' use by test statistic functions. 
#' 
#' @seealso Uses output from \link{formatClones} with similar arguments to 
#' \link{getTrees}. Output can be visualized with \link{plotTrees}, and tested
#' with \link{testPS}, \link{testSC}, and \link{testSP}.
#' 
#' @examples
#' \dontrun{
#' data(ExampleAirr)
#' ExampleAirr$sample_id <- sample(ExampleAirr$sample_id)
#' clones <- formatClones(ExampleAirr, trait="sample_id")
#' 
#' igphyml <- "~/apps/igphyml/src/igphyml"
#' btrees <- findSwitches(clones[1:2,], permutations=10, nproc=1,
#'    igphyml=igphyml, trait="sample_id")
#' plotTrees(btrees$trees[[4]])[[1]]
#' testPS(btrees$switches)
#' }
#' @export
findSwitches <- function(clones, permutations, trait, igphyml, 
                         fixtrees=FALSE, downsample=TRUE, tip_switch=20, nproc=1, 
                         dir=NULL, id=NULL, modelfile=NULL, build="pratchet", exec=NULL, 
                         quiet=0, rm_temp=TRUE, palette=NULL, resolve=2, rep=NULL,
                         keeptrees=FALSE, lfile=NULL, seq=NULL,
                         boot_part="locus", force_resolve=FALSE, ...){
  
  if(is.null(exec) && (!build %in% c("pratchet", "pml"))){
    stop("exec must be specified for this build option")
  }
  if(file.access(igphyml, mode=1) == -1) {
    stop("Igphyml executable at ", igphyml, " cannot be executed.")
  }
  if(downsample & !quiet & is.null(rep)){
    print(paste("Downsampling lineages to a maximum tip-to-switch ratio of",tip_switch))
  }
  if(!fixtrees & !quiet & is.null(rep)){
    print("Re-building bootstrapped trees. Use fixtrees=TRUE to use fixed topologies.")
    if(build=="igphyml"){
      stop("Bootstrapping while build=igphyml not yet supported (bootstraps sample nucleotide sites)")
    }
  }else if(is.null(rep)){
    print("Keeping tree topology constant. Use fixtrees=FALSE to bootstrap topologies.")
  }
  if(!is.null(dir)){
    dir <- path.expand(dir)
  }
  
  data <- clones$data
  if(!inherits(data, "list")){
    data <- list(data)
  }
  if(!inherits(data[[1]], "airrClone")){
    stop("Input data must be a list of airrClone objects")
  }
  if(fixtrees){
    if(!"trees" %in% names(clones)){
      stop("trees column must be included in input if fixtrees=TRUE (use getTrees first)")
    }
    if(!inherits(clones$trees[[1]], "phylo")){
      stop("Trees must be a list of class phylo")
    }
    trees <- clones$trees
  }else{
    trees <- NULL
  }
  if(is.null(id)){
    id <- "sample"
  }
  if(!is.null(dir)){
    if(!dir.exists(dir)){
      dir.create(dir)
    }
  }
  big <- FALSE
  if(sum(unlist(lapply(data, function(x)nrow(x@data)))) > 10000){
    big <- TRUE
  }
  if(!rm_temp && big){
    warning("Large dataset - best to set rm_temp=TRUE")
  }
  if(!is.null(igphyml)){
    igphyml <- path.expand(igphyml)
    if(!is.null(dir)){
      if(!dir.exists(dir)){
        dir.create(dir)
      }
    }else{
      dir <- alakazam::makeTempDir(id)
      if(big){
        warning("Large dataset - best to set dir and id params")
      }
    }
    if(is.null(trait)){
      stop("trait must be specified when running igphyml")
    }
    # remove problematic characters from trait values
    data <- lapply(data,function(x){
      x@data[[trait]] <- gsub(":|;|,|=| |_","-",x@data[[trait]])
      x})
        traits <- unique(unlist(lapply(data, function(x)unique(x@data[[trait]]))))
        if("N" %in% traits || "NTIP" %in% traits || "UCA" %in% traits){
            stop(paste("Forbidden trait values (e.g. N, NTIP, UCA) detected.",
                "Please remove or change, then re-run formatClones."))
        }

    if(is.null(modelfile)){
      states <- unique(unlist(lapply(data,function(x)x@data[,trait])))
      modelfile <- makeModelFile(states,file=file.path(dir,paste0(id,"_modelfile.txt")))
    }else{
      states <- readModelFile(modelfile)
    }
    if(sum(is.na(states) > 0)){
      stop("NA trait values detected, must be removed before trait analysis.")
    }
    # if trees are fixed, add trait value to tree tips
    if(fixtrees){
      for(i in 1:length(data)){
        da <- data[[i]]
        tr <- trees[[i]]
        tips <- tr$tip.label[tr$tip.label != "Germline"]
        m <- match(tips, da@data$sequence_id)
        if(sum(is.na(m)) > 0){
          stop(paste("Tips",tips[!tips %in% da@data$sequence_id],
                     "not found in clone object",da@clone))
        }
        tr$tip.label[tr$tip.label != "Germline"] <- 
          paste0(tips, "_", da@data[[trait]][m])
        trees[[i]] <- tr
      }
    }
    #if igphyml is specified, append trait value to sequence ids
    data <- lapply(data,function(x){
      x@data$sequence_id <- paste0(x@data$sequence_id,"_",x@data[[trait]])
      x})
  }
  if(build=="dnapars"){
    if(is.null(dir) || is.null(id)){
      stop("dir, and id parameters must be specified when running dnapars")
    }
  }
    if(is.null(rep)){
        reps <- as.list(1:permutations)
        l <- parallel::mclapply(reps,function(x)
            tryCatch(findSwitches(clones ,rep=x,
            trait=trait, modelfile=modelfile, build=build,
            exec=exec, igphyml=igphyml,
            id=id, dir=dir, permutations=permutations,
            nproc=1, rm_temp=rm_temp, quiet=quiet,
            fixtrees=fixtrees, resolve=resolve, keeptrees=keeptrees,
            lfile=lfile, seq=seq, downsample=downsample, tip_switch=tip_switch,
            boot_part=boot_part, force_resolve=force_resolve, ...),
            error=function(e)e),mc.cores=nproc)
        errors <- unlist(lapply(l, function(x) inherits(x, "error")))
        messages <- l[errors]
        if(sum(errors) > 0){
            me <- lapply(messages, function(x)warning(x$message))
            stop(paste("findSwitches failed for",sum(errors),
                "repetitions (see warnings, more info if nproc=1)"))
  }
    results <- list()
    results$switches <- NULL
    results$trees <- NULL
    if(!is.null(igphyml)){
      results$switches <- dplyr::bind_rows(lapply(l,function(x)x$switches))
    }
    if(keeptrees){
      results$trees <- lapply(l,function(x)x$trees)
    }
    if(rm_temp){
      if(file.exists(file.path(dir,paste0(id,"_modelfile.txt")))){
        unlink(file.path(dir,paste0(id,"_modelfile.txt")))
      }
    }
    return(results)
  }else{
    rm_dir=file.path(dir,paste0(id,"_recon_",rep))
    
    if(downsample){
      if(quiet > 3){print("downsampling clones")}
      rarefied <- lapply(1:length(data), function(x)
        downsampleClone(clone=data[[x]], 
                        tree=trees[[x]], trait=trait,
                        tip_switch=tip_switch))
      data <- lapply(rarefied, function(x)x$clone)
      seqs <- unlist(lapply(data, function(x)nrow(x@data)))
      index <- order(seqs, decreasing=TRUE)
      data <- data[index]
      if(fixtrees){
        trees <- lapply(rarefied, function(x)x$tree)
        trees <- trees[index]
      }
    }
    if(!fixtrees){
      temp_clones <- dplyr::tibble(data=data, clone_id = unlist(lapply(data, 
                            function(x)x@clone)), seqs = unlist(lapply(data,function(x)nrow(x@data))))
      #temp_clones <- getTrees(clones = temp_clones, build = build, nproc = nproc)
      clones_with_trees <- getBootstraps(clones = temp_clones, bootstraps = 1, nproc = nproc, 
                                         dir = dir, id = id, build = build, exec = exec,
                                         quiet = quiet, rm_temp = rm_temp, seq = seq,
                                         boot_part = boot_part, bootstrap_nodes = FALSE,
                                         switches = TRUE, ...)
      trees <- lapply(clones_with_trees$bootstrap_trees, function(x)x[[1]])
    }
    
    match <- unlist(lapply(1:length(data), function(x){
      data[[x]]@clone == trees[[x]]$name
    }))
    if(sum(!match) > 0){
      stop("Clone and tree names not in proper order!")
    }
    
    results <- list()
    if(!is.null(igphyml)){
      if(is.null(lfile)){    
        file <- writeLineageFile(data=data, trees=trees, dir=dir,
                                 id=id, trait=trait, rep=rep)
      }else{
        file <- lfile
      }
      rseed <- floor(stats::runif(1,0,10^7))+rep
      switches <- reconIgPhyML(file, modelfile, igphyml=igphyml, 
                               mode="switches", type="recon", id=rep, 
                               quiet=quiet, rm_files=FALSE, rm_dir=NULL, nproc=nproc,
                               resolve=resolve, rseed=rseed, force_resolve=force_resolve)
      permuted <- reconIgPhyML(file, modelfile, igphyml=igphyml, 
                               mode="switches", type="permute", id=rep, 
                               quiet=quiet, rm_files=FALSE, rm_dir=NULL, nproc=nproc,
                               resolve=resolve, rseed=rseed, force_resolve=force_resolve)
      if(!rm_temp){rm_dir=NULL}
      permuteAll <- reconIgPhyML(file, modelfile, igphyml=igphyml, 
                                 mode="switches", type="permuteAll", id=rep, 
                                 quiet=quiet, rm_files=rm_temp, rm_dir=rm_dir, nproc=nproc,
                                 resolve=resolve, rseed=rseed, force_resolve=force_resolve)
      switches <- rbind(switches,permuted)
      switches <- rbind(switches,permuteAll)
      switches$ID <- id
      #REP and CLONE are swapped up in the lower functions, this corrects them
      temp <- switches$REP
      switches$REP <- switches$CLONE
      switches$CLONE <- temp
      clone_index <- unlist(lapply(data,function(x)x@clone))
      switches$CLONE <- clone_index[switches$CLONE+1]
      results$switches <- switches
      
      # remove trait value from tips
      trees <- lapply(trees,function(x){
        ids <- strsplit(x$tip.label,split="_")
        x$tip.label <- unlist(lapply(ids,function(t)
          paste(t[1:(length(t)-1)],collapse="_")))
        x$tip.label[x$tip.label == x$name] <- "Germline"
        x
      })
    }
    if(keeptrees){
      results$trees <- trees
    }
    return(results)
  }
}

#' Deprecated! Please use findSwitches instead.
#' 
#' \code{bootstrapTrees} Phylogenetic bootstrap function.
#' @param clones         tibble \code{airrClone} objects, the output of 
#'                      \link{formatClones}
#' @param bootstraps    number of bootstrap replicates to perform
#' @param trait            trait to use for parsimony models (required if 
#'                      \code{igphyml} specified)
#' @param build           program to use for tree building (phangorn, dnapars)
#' @param exec           location of desired phylogenetic executable
#' @param igphyml        location of igphyml executible if trait models desired
#' @param id            unique identifer for this analysis (required if 
#'                      \code{igphyml} or \code{dnapars} specified)
#' @param dir           directory where temporary files will be placed (required
#'                      if \code{igphyml} or \code{dnapars} specified)
#' @param modelfile        file specifying parsimony model to use
#' @param fixtrees        keep tree topologies fixed?
#'                      (bootstrapping will not be perfomed)
#' @param nproc            number of cores to parallelize computations
#' @param quiet           amount of rubbish to print to console
#' @param rm_temp       remove temporary files (default=TRUE)
#' @param palette        deprecated
#' @param resolve        how should polytomies be resolved? 
#'                       0=none, 1=max parsminy, 2=max ambiguity + polytomy skipping,
#'                       3=max ambiguity
#' @param keeptrees     keep trees estimated from bootstrap replicates? (TRUE)
#' @param lfile         lineage file input to igphyml if desired (experimental)
#' @param rep             current bootstrap replicate (experimental)
#' @param seq           column name containing sequence information
#' @param downsample    downsample clones to have a maximum specified tip/switch ratio?
#' @param tip_switch    maximum allowed tip/switch ratio if downsample=TRUE
#' @param boot_part     is  "locus" bootstrap columns for each locus separately
#' @param force_resolve continue even if polytomy resolution fails?
#' @param ...        additional arguments to be passed to tree building program
#'
#' @return   A list of trees and/or switch counts for each bootstrap replicate.
#'  
#' @export
bootstrapTrees <- function(clones, bootstraps, nproc=1, trait=NULL, dir=NULL, 
                           id=NULL, modelfile=NULL, build="pratchet", exec=NULL, igphyml=NULL, 
                           fixtrees=FALSE,    quiet=0, rm_temp=TRUE, palette=NULL, resolve=2, rep=NULL,
                           keeptrees=TRUE, lfile=NULL, seq=NULL, downsample=FALSE, tip_switch=20,
                           boot_part="locus", force_resolve=FALSE,...){
  
  warning("boostrapTrees is depracated. Use findSwitches instead.")
  
  s = findSwitches(clones ,rep=rep, 
                   trait=trait, modelfile=modelfile, build=build, 
                   exec=exec, igphyml=igphyml, 
                   id=id, dir=dir, permutations=bootstraps,
                   nproc=1, rm_temp=rm_temp, quiet=quiet,
                   fixtrees=fixtrees, resolve=resolve, keeptrees=keeptrees,
                   lfile=lfile, seq=seq, downsample=downsample, tip_switch=tip_switch,
                   boot_part=boot_part, force_resolve=force_resolve, ...)
  return(s)
}


#' Get the tip labels as part of a clade defined by an internal node
#' 
#' \code{getSubTaxa} Gets the tip labels from a clade
#' @param  node    node number that defines the target clade
#' @param  tree    \code{phylo} object
#'
#' @return   A vector containing tip labels of the clade
#' @examples
#' # Get taxa from all subtrees
#' data(BiopsyTrees)
#' tree <- BiopsyTrees$trees[[8]]
#' all_subtrees <- lapply(1:length(tree$nodes), function(x)getSubTaxa(x, tree))
#' 
#' @export
getSubTaxa = function(node, tree){
  if(node > length(tree$tip.label) + tree$Nnode){
    stop(paste("Node", node, "too large for tree"))
  }
  children <- tree$edge[tree$edge[,1] == node, 2]
  if(length(children) == 0){
    if(node > length(tree$tip.label)){
      stop(paste("Malformed base case at node", node))
    }
    return(tree$tip.label[node])
  }
  tips <- unlist(lapply(children, function(x)getSubTaxa(x, tree)))
  return(tips)
}

# KEN: Try to keep lines to 80 characters or less, which is here -------------->|

# Turn your tree data into a nodes based dataframe
# 
# \code{lones} Filler
# @param    input_tree     \code{phylo} object
# @param    bootstrap_number      Which bootstrap replicate to use. With only 
# @param                          one tree, this has to be one. 
#
# @return   A dataframe that lists out the value of found or absent tips for each 
#           node within a tree. This is done for each tree inputted.
# bootstrapped sequences
splits_func <- function(input_tree, bootstrap_number){
  tree <- input_tree[[bootstrap_number]] 
  splits <- data.frame(found=I(lapply(
    (length(tree$tip.label) + 1):(length(tree$tip.label) +  
                                    dplyr::n_distinct(tree$edge[,1])),
    function(x)getSubTaxa(x, tree))))
  splits$node <- (length(tree$tip.label) + 1):(length(tree$tip.label) +
                                              dplyr::n_distinct(tree$edge[,1]))
  # find the difference between tip labels and the tips in 'found'
  full_tips <- tree$tip.label 
  absent <- (lapply(1:length(splits$found), 
                    function(x)dplyr::setdiff(full_tips, splits$found[[x]])))
  splits <- cbind(splits, data.frame(absent=I(absent)))
  splits$tree_num <- bootstrap_number
  # reorder it to make sense -- tree number, node number, found tips, and absent
  splits <- splits[, c(4, 2, 1, 3)]
  return(splits)
}

# Match the tips found in the various nodes of two different trees. 
# 
# \code{lones} Filler
# @param    tree_comp_df     The tree the bootstrap tree nodes should be 
#                            compared to. This needs to already be a dataframe 
#                            made from splits_func.
# @param    bootstrap_df     The dataframe made by the splits_func that has the 
#                            found and absent tips in the bootstrap trees.
# @param    nproc            The number of processors to be used
#
# @return   Returns a vector with the number of matches. 
matching_function_parallel <- function(tree_comp_df, bootstrap_df, nproc){
  match_vector <- c()
  match_vector = parallel::mclapply(unique(tree_comp_df$node), function(node){
    # KEN: There must be a more efficient way of doing this but I can't think of 
    # one right now
    sub_full_tree_df <- tree_comp_df[tree_comp_df$node==node,]
    matches = unlist(lapply(1:length(bootstrap_df$tree_num), function(x){
      match_test1 <- setequal(bootstrap_df$found[[x]], sub_full_tree_df$found[[1]])
      match_test2 <- setequal(bootstrap_df$found[[x]], sub_full_tree_df$absent[[1]])
      if(match_test1 || match_test2){
        return(1)
      }else{
        return(0)
      }}))
    dplyr::tibble(nodes = node, matches = sum(matches))
  },mc.cores=nproc)
  match_vector <- dplyr::bind_rows(match_vector)
  return(match_vector)
}


# KEN: Output very suspicious, this may be one we need to do ourselves..
# Build a bootstrap consensus tree using list of bootstrapped trees
# 
# \code{lones} Filler
# @param    trees     List of trees to use for tree building
#
# @return   Returns a consensus tree. 
consensus_tree <- function(trees){
  consensus <- ape::consensus(trees, check.labels=TRUE)
  consensus$edge.length <- rep(1, nrow(consensus$edge))
  consensus <- rerootTree(consensus, germline="Germline", verbose=0)
  return(consensus)
}


# Build some trees based on the clones object's $data and other parameters.
# This results in a list of trees. 
#
# \code{lones} Filler
# @param data          an AIRRE clone's data object (e.g. clones$data)
# @param seq           column name containing sequence information
# @param build         program to use for tree building (phangorn, dnapars, igphyml)
# @param boot_part     is  "locus" bootstrap columns for each locus separately
# @param exec          location of desired phylogenetic executable
# @param dir           directory where temporary files will be placed (required
#                      if \code{igphyml} or \code{dnapars} specified)
# @param rm_temp       remove temporary files (default=TRUE)
# @param rm_dir        remove temporary directory (default=TRUE)
# @param id            unique identifier for this analysis (required if 
#                      \code{igphyml} or \code{dnapars} specified)
# @param quiet           amount of rubbish to print to console
# @param rep             current bootstrap replicate (experimental)
# @param by_codon      Bootstrap by codon 
# @param starting_tree A starting tree 
# @return              Returns a list of trees.
makeTrees <- function(clones, seq, build, boot_part, exec, dir, rm_temp=TRUE, id, 
                      quiet=0, rep=1, by_codon = TRUE, starting_tree=NULL, ...){
  if(quiet > 0){
    print(paste0("Making trees for rep ", rep))
  }
  #KBH in future iterations we'll want to make the boostrapping optional, but it's okay for now
  data <- clones$data
  # randomize the data_tmp file -- shuffle the rows CGJ 4/13/23
  for(i in 1:length(data)){
    data[[i]]@data <- data[[i]]@data[sample(1:nrow(data[[i]]@data), replace = FALSE),]
  }
  data_tmp <- list()
  for(i in 1:length(data)){
    data_tmp[[i]] <- bootstrapClones(data[[i]], reps = 1, 
                                     partition = boot_part, by_codon = by_codon)[[1]]
  }
  reps <- as.list(1:length(data_tmp))
  if(is.null(seq)){
    seqs <- unlist(lapply(data_tmp,function(x)x@phylo_seq))
  }else{
    seqs <- rep(seq,length=length(data_tmp))
  }
  if(!is.null(starting_tree)){
    from_getTrees <- TRUE
  } else{
    from_getTrees <- FALSE
  }
  if(build=="pratchet"){
    trees <- lapply(reps,function(x)
      tryCatch(buildPratchet(data_tmp[[x]],seq=seqs[x],...),error=function(e)e))
    trees <- list(trees)
  } else if(build=="dnapars" || build=="dnaml"){
    trees <- lapply(reps,function(x)
      tryCatch(buildPhylo(data_tmp[[x]],
                          exec=exec,
                          temp_path = file.path(dir,paste0(id,"_", rep, "_trees_",x)),
                          rm_temp = rm_temp,
                          seq=seqs[x], ...), error=function(e)e))
    trees <- list(trees)
  }else if(build=="pml"){
    trees <- lapply(reps,function(x)tryCatch(buildPML(data_tmp[[x]],seq=seqs[x],
                                quiet=quiet, rep=rep,...), error=function(e)e))
    
    trees <- list(trees)
  } else if(build=="raxml"){ # CGJ 2/20/23
    trees <- lapply(reps,function(x)tryCatch(buildRAxML(data_tmp[[x]],
                          seq=seqs[x],
                          exec = exec,
                          trees = starting_tree[[x]],
                          rm_files = rm_temp,
                          rep = rep,
                          dir = dir,
                          rseed = x,
                          from_getTrees=from_getTrees,...),error=function(e)e))
  } else if(build=="igphyml"){
    if(rm_temp){
      rm_dir <- file.path(dir,paste0(id,rep))
    }else{
      rm_dir <- NULL
    }
    trees <- 
      tryCatch(buildIgphyml(data_tmp, 
                            igphyml = exec,
                            temp_path = file.path(dir,paste0(id,rep)),
                            rm_files=rm_temp,
                            rm_dir = rm_dir,
                            id=id, ...), error=function(e)e)
    trees <- list(trees)
    if(inherits(trees, "error")){
      stop(trees)
    }
    
  }else{
    stop("build specification ",build," not recognized")
  } 
  
  if(build != "igphyml"){
    errors <- unlist(lapply(trees, function(x) inherits(x, "error")))
    messages <- trees[errors]
    errorclones <- clones$clone_id[errors]
    trees <- trees[!errors]
    if(length(errorclones) > 0){
      warning(paste("Tree building failed for clones",
                    paste(errorclones,collapse=", ")))
      me <- lapply(messages, function(x)warning(x$message))
    }
  }
  return(trees)
}


#' Creates a bootstrap distribution for clone sequence alignments, and returns  
#' estimated trees for each bootstrap replicate as a nested list as a new input 
#' tibble column.
#' 
#' \code{getBootstraps} Phylogenetic bootstrap function.
#' @param clones         tibble \code{airrClone} objects, the output of 
#'                      \link{formatClones}
#' @param bootstraps    number of bootstrap replicates to perform
#' @param nproc            number of cores to parallelize computations
#' @param bootstrap_nodes a logical if the the nodes for each tree in the trees
#'                        column (required) should report their bootstrap value
#' @param dir           directory where temporary files will be placed (required
#'                      if \code{igphyml} or \code{dnapars} specified)
#' @param id            unique identifer for this analysis (required if 
#'                      \code{igphyml} or \code{dnapars} specified)
#' @param build           program to use for tree building (phangorn, dnapars, igphyml)
#' @param exec           location of desired phylogenetic executable
#' @param quiet           amount of rubbish to print to console
#' @param rm_temp       remove temporary files (default=TRUE)
#' @param rep             current bootstrap replicate (experimental)
#' @param seq           column name containing sequence information
#' @param boot_part     is  "locus" bootstrap columns for each locus separately
#' @param by_codon      a logical if the user wants to bootstrap by codon or by 
#'                      nucleotide. Default (codon based bootstrapping) is TRUE.
#' @param starting_tree An indicator to use the existing trees column as the starting trees for RAxML
#' @param switches   a logical indicator to allow findSwitches to do permutations. 
#' @param ...        additional arguments to be passed to tree building program
#'
#' @return   The input clones tibble with an additional column for the bootstrap replicate trees.
#'  
#' @export
getBootstraps <- function(clones, bootstraps,
                          nproc=1, bootstrap_nodes=TRUE, dir=NULL, id=NULL, build="pratchet", 
                          exec=NULL, quiet=0, rm_temp=TRUE, rep=NULL, seq=NULL,
                          boot_part="locus", by_codon = TRUE, starting_tree=FALSE,
                          switches=FALSE, ...){
  if(is.null(exec) && (!build %in% c("pratchet", "pml"))){
    stop("exec must be specified for this build option")
  }
  if(build=="igphyml" && file.access(exec, mode=1) == -1) {
    stop("Igphyml executable at ", exec, " cannot be executed.")
  }
  if(!is.null(dir)){
    dir <- path.expand(dir)
  }
  
  data <- clones$data
  if(is.null(id)){
    id <- "sample"
  }
  if(!inherits(data,"list")){
    data <- list(data)
  }
  if(!is.null(dir)){
    if(!dir.exists(dir)){
      dir.create(dir)
    }
  }
  if(!inherits(data[[1]], "airrClone")){
    stop("Input data must be a list of airrClone objects")
  }
  big <- FALSE
  if(sum(unlist(lapply(data, function(x)nrow(x@data)))) > 10000){
    big <- TRUE
  }
  if(!rm_temp && big){
    warning("Large dataset - best to set rm_temp=TRUE")
  }
  if(build == "igphyml"){
    # CGJ 9/5/23 checkUsagePackage("dowser") found it
    # igphyml <- path.expand(exec)
    if(!is.null(dir)){
      if(!dir.exists(dir)){
        dir.create(dir)
      }
    } else{
      dir <- alakazam::makeTempDir(id)
      if(big){
        warning("Large dataset - best to set dir and id params")
      }
    }
  }
  if(build=="dnapars" || build=="dnaml"){
    exec <- path.expand(exec)
    if(!is.null(dir)){
      if(!dir.exists(dir)){
        dir.create(dir)
      }
    } else{
      dir <- alakazam::makeTempDir(id)
      if(big){
        warning("Large dataset - best to set dir and id params")
      }
    }
  }
  if(!switches){
    if(bootstrap_nodes){
      if(!"trees" %in% colnames(clones)){
        stop("A trees column created by using getTrees() is required if 
           bootstrap_nodes=TRUE")
      }
    }
    build_used <- gsub("phangorn::", "", clones$trees[[1]]$tree_method)
    build_used <- gsub("phylip::", "", build_used)
    build_used <- gsub("\\:.*", "", build_used)
    build_used <- gsub("optim.", "", build_used)
    if(grepl("igphyml", build_used)){
      build_used <- "igphyml"
    }
    if(grepl("RAxML", build_used)){
      build_used <- "raxml"
    }
    if(build_used == "prachet"){
      build_used <- "pratchet"
    }
    if(build != build_used){
      stop(paste0("Trees and bootstrapped trees need to be made using the same method.",
                  " Use the same build option in getTrees as getBootstraps.",
                  " getBootstraps is trying to use a ", build, 
                  " build, but getTrees used ", build_used, " to build trees."))
    }
    # CGJ 2/20/23
    if(starting_tree){
      if("trees" %in% colnames(clones)){
        starting_tree <- clones$trees
      } else{
        stop("starting_trees cannot be set as TRUE without an already made trees",
             " column. Use getTrees() to get the trees column.")
      }
    } else{
      starting_tree <- NULL
    }
  }
  if(build != "igphyml" | build != "raxml"){
    bootstrap_trees <- unlist(parallel::mclapply(1:bootstraps, function(x)
      tryCatch(makeTrees(clones=clones, seq=seq, build=build, boot_part=boot_part,
                         exec=exec, dir=dir, rm_temp=rm_temp, id=id, quiet=quiet, 
                         rep=x, by_codon = by_codon), 
               error=function(e)e), mc.cores=nproc), recursive = FALSE)
  }else if(build == "igphyml" | build == "raxml"){
    bootstrap_trees <- unlist(parallel::mclapply(1:bootstraps, function(x)
      tryCatch(makeTrees(clones=clones, seq=seq[x], build=build, boot_part=boot_part,
                         exec=exec, dir=dir, rm_temp=rm_temp, id=id, quiet=quiet, 
                         rep=x, by_codon = by_codon, starting_tree = starting_tree), 
               error=function(e)e), mc.cores=nproc), recursive = FALSE)
  }
  clones$bootstrap_trees <- lapply(1:nrow(clones), function(x)list())
  # make raxml format friendly CGJ 2/20/23
  if(build == "raxml"){
    new_boots <- c()
    for(i in 1:bootstraps){
      values <- 1:nrow(clones)*i
      if(i == 1){
        tobind <- bootstrap_trees[min(values):max(values)]
      } else{
        tobind <- bootstrap_trees[((max(values)/min(values))+1):max(values)]
      }
      new_boots <- append(new_boots, list(tobind))
    }
    bootstrap_trees <- new_boots
  }
  for(i in 1:length(clones$clone_id)){
    clones$bootstrap_trees[[i]] <- lapply(bootstrap_trees, function(x)x[[i]])
  }
  errors <- c()
  messages <- c()
  for(clone in unique(clones$clone_id)){
    if(quiet > 0){
      print(clone)
    }
    sub_clone <- dplyr::filter(clones, !!rlang::sym("clone_id") == clone)
    tree_error <- c()
    tree_message <- c()
    for(i in 1:bootstraps){
      if(inherits(sub_clone$bootstrap_trees[[1]][[i]], "error")){
        tree_error <- append(tree_error, clone)
        tree_message <- append(tree_message, sub_clone$bootstrap_trees[[1]][[i]]$message)
      }
    }
    errors <- append(errors, tree_error)
    messages <- append(messages, tree_message)
  }
  
  for(clone in unique(errors)){
    idx <- which(errors == clone)
    wow <- messages[idx]
    warning(paste0("Clone ", clone, " has been REMOVED from the clones object ",
                   "because it failed to properly bootstrap due in at least one iteration of ",
                   "bootstrapping due to ", unique(wow), "."))
  }
  clones <- dplyr::filter(clones, !(!!rlang::sym("clone_id") %in% unique(errors)))
  if(!bootstrap_nodes){
    return(clones)
  } else{
    for(clone in 1:length(clones$clone_id)) {
      if(quiet >0){
        print(paste0("bootstrapping clone ", clones$clone_id[clone]))
      }
      b_trees <- clones$bootstrap_trees[[clone]] # KEN: safer if clones is empty
      tree_comp_df <- splits_func(list(clones$trees[[clone]]), 1)
      bootstraps_df <- lapply(1:bootstraps, function(x)splits_func(b_trees,x))
      bootstraps_df <- do.call(rbind, bootstraps_df)
      matches_df <- matching_function_parallel(tree_comp_df, bootstraps_df, nproc)
      for(node in 1:max(matches_df$nodes)){
        if(node %in% setdiff(clones$trees[[clone]]$edge[,2], clones$trees[[clone]]$edge[,1])){
          clones$trees[[clone]]$nodes[[node]]$bootstrap_value <- NA
        } else{
          clones$trees[[clone]]$nodes[[node]]$bootstrap_value <- 
            dplyr::filter(matches_df, !!rlang::sym("nodes") == node)$matches
        }
      }
      if(quiet >0){
        print(paste0("Clone ", clones$clone_id[clone], " bootstrapping completed"))
      }
    }
  }
  return(clones)
}


#' Exports the phylogentic trees from the airrClone object
#' 
#' \code{exportTrees}   Exports phylogenetic trees
#' @param clones         tibble \code{airrClone} objects, the output of 
#'                      \link{formatClones}
#' @param file       The file path and name of where the trees will be saved
#' @param ...        additional arguments to be passed
#'  
#' @export
exportTrees <- function(clones, file, tree_column = "trees", ...){
  # check to see if the trees column is there 
  if(alakazam::checkColumns(clones, tree_column)){
    ape::write.tree(phy = clones$trees, file = file, ...)
  } else{
    stop(paste(tree_column, "not found in the input airrClone object. Please",
               "specify what column cotains the phylogentic trees."))
  }
  
}

#' Write the sequences used in tree building to a fasta format. If there are more 
#' than one tree in airrClone output the sequence id will be followed by "|clone_id".
#' 
#' \code{writeCloneSequences}   Exports the sequences used in tree building. 
#' @param clones         tibble \code{airrClone} objects, the output of 
#'                      \link{formatClones}
#' @param file       The file path and name of where the sequences will be saved
#'  
#' @export
writeCloneSequences <- function(clones, file){
  for(i in 1:nrow(clones)){
    clone_id <- clones$clone_id[i]
    # grab the germline 
    if(clones$data[[i]]@phylo_seq == "sequence"){
      germline <- clones$data[[i]]@germline
    } else if(clones$data[[i]]@phylo_seq == "hlsequence"){
      germline <- clones$data[[i]]@hlgermline
    } else if(clones$data[[i]]@phylo_seq == "lsequence"){
      germline <- clones$data[[i]]@lgermline
    }
    write(paste0(">Germline|", clone_id), file, append = TRUE)
    write(germline, file, append = TRUE)
    for(j in 1:nrow(clones$data[[i]]@data)){
      seq_id <- clones$data[[i]]@data$sequence_id[j]
      if(clones$data[[i]]@phylo_seq == "sequence"){
        sequence <- clones$data[[i]]@data$sequence[j]
      } else if(clones$data[[i]]@phylo_seq == "hlsequence"){
        sequence <- clones$data[[i]]@data$hlsequence[j]
      } else if(clones$data[[i]]@phylo_seq == "lsequence"){
        sequence <- clones$data[[i]]@data$lsequence[j]
      }
      write(paste0(">", seq_id, "|", clone_id), file, append = TRUE)
      write(sequence, file, append = TRUE)
    }
  }
}