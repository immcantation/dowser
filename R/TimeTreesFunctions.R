# Functions for performing time tree analysis on B cell lineage trees

#' Iteratively resume getTimeTrees until convergence, as defined by 
#' all parameters (except those in \code{ignore} vector) having ESS 
#' greater than or equal to the specified ess_cutoff
#'
#' \code{getTimeTreesIterate} Iteratively resume getTimeTrees til convergence.
#' @param    clones     a tibble of \code{airrClone} objects, the output of
#'                      \link{formatClones}
#' @param    iterations Maximum number of iterations
#' @param    ess_cutoff Minimum number of ESS for all parameters
#' @param    ignore     Vector of parameters to ignore for ESS calculation
#' @param    quiet      quiet notifications if > 0
#' @param    ...        Additional arguments for getTimeTrees
#'
#' @return   A tibble of \code{tidytree} and \code{airrClone} objects.
#'
#' @details
#' For examples and vignettes, see https://dowser.readthedocs.io
#'
#' @export
getTimeTreesIterate <- function(clones, iterations=10, ess_cutoff=200,
  ignore = c("traitfrequencies"), quiet=0, ...){

  resume = NULL
  iter = 0
  while((length(resume) > 0 || iter == 0 ) && iter < iterations){
        
        if(quiet < 1){
          print(paste("Starting iteration", iter))
        }

        clones = getTimeTrees(clones, resume_clones=resume, ...)

        params = clones$parameters
        for(regex in ignore){
            params = lapply(params, function(x){
                dplyr::filter(x, !grepl(regex, !!rlang::sym("item")))
                })
        }
        
        clones$below_ESS = sapply(params, 
          function(x)sum(x$ESS[!x$item %in% ignore] < ess_cutoff, na.rm=TRUE))
        if(quiet < 1){
          print(clones$below_ESS)
          ess_items = unlist(sapply(params, function(x)x$item[x$ESS[!x$item %in% ignore] < ess_cutoff]))
          if(length(ess_items) > 0){
            print(table(ess_items))
          }
        }
    
        resume = dplyr::filter(clones, !!rlang::sym("below_ESS") > 0)$clone_id
        iter = iter + 1
    }
    if(iter == iterations & length(resume) != 0){
      warning(paste(paste(resume, collapse=","), "failed to converge after",
        iterations, "iterations"))
    }
    return(clones)
}

#' Estimate time trees by running BEAST on each clone
#' Applies XML \code{template} to each clone
#'
#' \code{getTimeTrees} Tree building function.
#' @param    clones     a tibble of \code{airrClone} objects, the output of
#'                      \link{formatClones}
#' @param    template   XML template
#' @param    beast      location of beast binary directory (beast/bin)
#' @param    dir        directory where temporary files will be placed.
#' @param    id         unique identifer for this analysis
#' @param    mcmc_length  Number of MCMC iterations
#' @param    time         Name of sample time column  
#' @param    log_every    Frequency of states logged. "auto" will divide
#'                        mcmc_length by log_target         
#' @param    burnin       Burnin percent (default 10)                 
#' @param    trait        Trait coolumn used         
#' @param    resume_clones  Clones to resume for mcmc_length more iterations            
#' @param    include_germline Include germline in analysis?     
#' @param    seq          Sequence column in data      
#' @param    germline_range   Possible date range of germline tip     
#' @param    java         Use the -java flag for BEAST run      
#' @param    seed         Used for the -seed option for BEASTrun    
#' @param    log_target   Target number of samples over mcmc_length         
#' @param    tree_states  Use \code{states} vector for starting tree
#' @param    nproc      Number of cores for parallelization. Uses 1 core per tree.
#' @param    quiet      amount of rubbish to print to console
#' @param    rm_temp    remove temporary files (default=TRUE)
#' @param    trees      optional list of starting trees, either phylo objects or newick strings
#' @param    ...        Additional arguments passed to tree building programs
#'
#' @return   A list of \code{phylo} objects in the same order as \code{data}.
#'
#' @details
#' For examples and vignettes, see https://dowser.readthedocs.io
#'
#' @seealso \link{getTrees}, \link{readBEAST}
#' @export
getTimeTrees <- function(clones, template, beast, dir, id, time,  
        mcmc_length=30000000, log_every="auto", 
        burnin=10, trait=NULL, resume_clones=NULL, nproc=1, quiet=0, 
        rm_temp=FALSE, include_germline=TRUE, seq="sequence", 
        germline_range=c(-10000,10000), java=TRUE, seed=NULL, log_target=10000, 
        tree_states=FALSE, trees=NULL, ...){

  if(is.null(beast)){
    stop("BEAST bin directory must be specified for this build option")
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
  if(sum(clones$clone_id != sapply(clones$data, function(x)x@clone)) > 0){
    stop("clone_id and airrClone values not identical")
  }
  if(!is.null(resume_clones) && !"trees" %in% names(clones)){
    stop("trees column not found. resume_clones can only be used on data after at least one getTimeTrees run.")
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
  if(is.null(id)){
    id <- "sample"
  }
  if(sum(unlist(lapply(data, function(x)nrow(x@data) > 100))) > 0 & quiet < 1){
    warning("Some clones contain > 100 sequences, may be slow")
  }

  if(!time %in% names(data[[1]]@data)){
    stop(paste(time, "column not found in data"))
  }
  if(!"numeric" %in% class(data[[1]]@data[[time]])){
    stop("time column must be numeric")
  }
  if(is.null(dir)){
    stop("dir must be specified when running BEAST")
  }
#  if(!is.null(time)){
#    # remove problematic characters from trait values
#    capture <- lapply(data,function(x)
#      if(sum(is.na(as.numeric(x@data[[time]]))) > 0){
#        stop(paste("clone",x@clone,
#          ": trait values must be non-missing numeric values when running beast"))
#      })
#  }
  if(!is.null(dir)){
    if(!dir.exists(dir)){
      dir.create(dir)
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
    rm_dir <- file.path(dir,paste0(id,"_beast"))
  }

  reps <- as.list(1:length(data))
  if(is.null(seq)){
    seqs <- unlist(lapply(data,function(x)x@phylo_seq))
  }else{
    seqs <- rep(seq,length=length(data))
  }
  
  # write buildBeast
  trees <- tryCatch(buildBeast(data,
                            time=time,
                            trait=trait,
                            mcmc_length=mcmc_length,
                            burnin=burnin,
                            beast=beast,
                            dir=dir,
                            id=id,
                            nproc=nproc,
                            template=template,
                            include_germline=include_germline,
                            resume_clones=resume_clones, 
                            log_every=log_every,
                            germline_range=germline_range,
                            java=java,
                            seed=seed,
                            log_target=log_target,
                            tree_states=tree_states,
                            trees=trees,
                            ...
                            ),error=function(e)e)

  if(inherits(trees, "error")){
    stop(trees)
  }
  
  if(length(trees) == 0){
    stop("No trees left!")
  }

  if(is.null(resume_clones)){
    # make sure trees, data, and clone objects are in same order
    tree_names <- unlist(lapply(trees, function(x)x@info$name))
    data_names <- unlist(lapply(data, function(x)x@clone))
    m <- match(tree_names, data_names)
    data <- data[m]
  
    # Sanity check
    match <- unlist(lapply(1:length(data), function(x){
      data[[x]]@clone == trees[[x]]@info$name
    }))
    if(sum(!match) > 0){
      stop("Clone and tree names not in proper order!")
    }
    m <- match(tree_names, clones$clone_id)
    clones <- clones[m,]
    clones$trees <- trees
  }else{
    tree_names <- unlist(lapply(trees, function(x)x@info$name))
    m <- match(tree_names, clones$clone_id)
    clones$trees[m] <- trees
  }

  #Sanity check
  if(sum(clones$clone_id != 
         unlist(lapply(clones$trees,function(x)x@info$name))) > 0){
    stop("Tree column names don't match clone IDs")
  }
  clones$parameters <- lapply(clones$trees, function(x)x@info$parameters)
  clones$ESS_100 = sapply(clones$parameters, function(x)sum(x$ESS < 100, na.rm=TRUE))
  clones$ESS_200 = sapply(clones$parameters, function(x)sum(x$ESS < 200, na.rm=TRUE))

  clones
}


#' Read in a directory from a BEAST run. Runs treeannotator and loganalyser.
#' 
#' @param    data     a list of \code{airrClone} objects
#' @param    template   XML template
#' @param    beast      location of beast binary directory (beast/bin)
#' @param    dir        directory where temporary files will be placed.
#' @param    id         unique identifer for this analysis
#' @param    mcmc_length  Number of MCMC iterations
#' @param    time         Name of sample time column 
#' @param    log_every    Frequency of states logged. \code{auto} will divide mcmc_length by log_target         
#' @param    burnin       Burnin percent (default 10)                 
#' @param    trait        Trait coolumn used         
#' @param    asr          Log ancestral sequences?
#' @param    full_posterior  Read un full distribution of parameters and trees?
#' @param    resume_clones  Clones to resume for mcmc_length more iterations            
#' @param    include_germline Include germline in analysis?     
#' @param    start_date       Starting date of time tree if desired
#' @param    max_start_date   Maximum starting date of time tree if desired
#' @param    germline_range   Possible date range of germline tip     
#' @param    java         Use the -java flag for BEAST run      
#' @param    seed         Used for the -seed option for BEASTrun    
#' @param    log_target   Target number of samples over mcmc_length         
#' @param    tree_states  Use \code{states} vector for starting tree
#' @param    nproc      Number of cores for parallelization. Uses 1 core per tree.
#' @param    quiet      amount of rubbish to print to console
#' @param    low_ram    run with less memory (slower)  
#' @param    trees                    optional list of starting trees, either phylo objects or newick strings
#' @param    start_edge_length        edge length to use for all branches in starting tree 
#' @param    ...      Additional arguments for XML writing functions
#'
#' @return   The input clones tibble with an additional column for the bootstrap replicate trees.
#'  
#' @seealso \link{getTimeTrees}
#' @export
buildBeast <- function(data, beast, time, template, dir, id, mcmc_length = 1000000, 
                   resume_clones=NULL, trait=NULL, asr=FALSE,full_posterior=FALSE,
                   log_every="auto",include_germline = TRUE, nproc = 1, quiet=0, 
                   burnin=10, low_ram=TRUE, germline_range=c(-10000,10000), java=TRUE, 
                   seed=NULL, log_target=10000, trees=NULL, tree_states=FALSE, 
                   start_edge_length=100, start_date=NULL, max_start_date=NULL,...) {

  beast <- path.expand(beast)
  beast_exec <- file.path(beast,"beast")
  if(file.access(beast_exec, mode=1) == -1) {
    stop("The file ", beast_exec, " cannot be executed.")
  }
  annotator_exec <- file.path(beast,"treeannotator")
  if(file.access(annotator_exec, mode=1) == -1) {
    stop("The file ", annotator_exec, " cannot be executed.")
  }
  analyser_exec <- file.path(beast,"loganalyser")
  if(file.access(analyser_exec, mode=1) == -1) {
    stop("The file ", analyser_exec, " cannot be executed.")
  }
  if(burnin > 100 || burnin < 0){
    stop("burnin must be between 0 and 100 (represents %)")
  }

  # setting log_every to get log_target samples given chain length
  if(log_every == "auto"){
    log_every <- max(floor(mcmc_length/log_target), 1)
  }
  if(!grepl("2\\.7", beast_exec) && quiet < 1){
    warning("most templates only compatible with only beast 2.7")
  }

  trait_list <- NULL
  if(!is.null(trait)){
    trait_list <- sort(unique(unlist(lapply(data, function(x)x@data[[trait]]))))
  }

  if(!is.null(resume_clones)){
    not_found = resume_clones[!resume_clones %in% sapply(data, function(x)x@clone)]
    if(length(not_found) > 0){
      stop("Clones ", paste(not_found, collapse=","), " not in data")
    }
    cat("Re-running clones", paste(resume_clones, collapse=","), "\n")
    data <- data[sapply(data, function(x)x@clone %in% resume_clones)]
    if(!is.null(trees)){
      trees <- trees[sapply(trees, function(x)x$name %in% resume_clones)]
    }
  }
  if(!is.null(trees)){
    check <- sapply(1:length(data), function(x){
      data[[x]]@clone == trees[[x]]$name
    })
    if(sum(!check) > 0){
      stop("Clone and starting tree names do not match")
    }
  }
  
  # Create the XML file using the current template and include germline
  xml_filepath <- write_clones_to_xmls(data, 
      mcmc_length=mcmc_length,
      log_every=log_every,
      id=id,
      outfile=file.path(dir, id), 
      time=time, 
      trait=trait, 
      trait_list=trait_list,
      template=template,
      include_germline_as_root=include_germline,
      include_germline_as_tip=include_germline, 
      germline_range=germline_range,
      tree_states=tree_states, 
      start_edge_length=start_edge_length,
      trees=trees,
      start_date=start_date, 
      max_start_date=max_start_date,
      ...)

  xml_filepath <- xml_filepath[!is.na(xml_filepath)]

  # Run BEAST on each tree sequentially
  # TODO: option to parallelize by tree?
  capture <- parallel::mclapply(1:length(xml_filepath), function(x) {
    y <- xml_filepath[x]
    overwrite <- "-overwrite"
    if(!is.null(resume_clones)){
      overwrite <- "-resume"
    }
    if(java){
      command <- paste0(
      "\ ", "-threads\ ", 1,
      "\ ", "-working\ ", 
      "\ ", "-java\ ", 
      "\ ",overwrite, "\ ")
    }else{
      command <- paste0(
      "\ ", "-threads\ ", 1,
      "\ ", "-working\ ", 
      "\ ",overwrite, "\ ")
    }

    if(is.null(seed)){
      command <- paste0(command, y)
    }else{
      command <- paste0(command, "-seed ",seed, "\ ",y)
    }

    console_out <- paste(gsub(".xml$","_console.log",y))
    
    if(is.null(resume_clones)){  
      command <- paste0(command, " > ", console_out)
    }else{
      command <- paste0(command, " >> ", console_out)
    }
    
    if(quiet < 1){
      print(paste(beast_exec,command))
    }
      
    params <- list(beast_exec, command, stdout=TRUE, stderr=TRUE)

    status <- tryCatch(do.call(base::system2, params), error=function(e){
         print(paste("BEAST error: ",e));
         return(e)
     }, warning=function(w){
         print(paste("BEAST warnings ",w));
         return(w)
     })
    status
    }, mc.cores=nproc)

  for(i in 1:length(capture)){
    if("error" %in% class(capture[[i]])){
      print(capture[[i]])
      stop(paste("Error running BEAST (see above), clone", data[[i]]@clone))
    }
  }

 trees <- readBEAST(clones=data, dir=dir, id=id, beast=beast, burnin=burnin, 
  trait=trait, quiet=quiet, nproc=nproc, full_posterior=full_posterior, asr=asr, 
  low_ram=low_ram)

  return(trees)
}


#' Takes an airr clone object and returns BEAST2 Alignment xml of the sequences
#' 
#' @param    clone                    an \code{airrClone} object
#' @param    id                       unique identifer for this analysis
#' @param    include_germline_as_tip  include the germline as a tip in the alignment?
#'
#' @return   String of BEAST2 Alignment and TaxonSet xml
#'  
create_alignment <- function(clone, id, include_germline_as_tip) {
  
  all_seqs <- ""
  
  for (i in 1:nrow(clone@data)) {
    # create a sequence object
    sequence <- clone@data[i, ]
    sequence_xml <- 
      paste0('\t<sequence id="seq_', sequence$sequence_id, 
             '" spec="Sequence" taxon="', sequence$sequence_id, 
             '" totalcount="4" value="', sequence$sequence, '" />\n')
    all_seqs <- paste0(all_seqs, sequence_xml)
  }
  
  if (include_germline_as_tip) {
    germline_sequence_xml <- 
      paste0('\t<sequence id="seq_', 'Germline', 
             '" spec="Sequence" taxon="', 'Germline', 
             '" totalcount="4" value="', clone@germline, '" />\n')
    all_seqs <- paste0(all_seqs, germline_sequence_xml)
  }
  
  alignment_xml <- 
    paste0('<data id="', id, "_", clone@clone, '" spec="Alignment" name="alignment">\n', 
           all_seqs, 
           '</data>')
  
  # create a taxon set for this alignment
  taxon_set <- 
    paste0('<taxa id="TaxonSet.', id, "_", clone@clone, '" spec="TaxonSet">\n', 
           '\t<alignment idref="', id, "_", clone@clone, '"/>\n', 
           '</taxa>')
  
  return(paste0(alignment_xml, '\n', taxon_set, '\n'))
}


#' Takes an airr clone object and returns BEAST2 rootfreqs xml of the germline
#' 
#' @param    clone                    an \code{airrClone} object
#' @param    id                       unique identifer for this analysis
#'
#' @return   String of XML setting the root frequencies to the germline sequence
#'  
create_root_freqs <- function(clone, id) {
  if (any(grepl("N", clone@germline))) {
    freqs <- clone@germline
    freqs <- gsub("N", " 0.25,0.25,0.25,0.25;", freqs)
    freqs <- gsub("A", " 1,0,0,0;", freqs)
    freqs <- gsub("C", " 0,1,0,0;", freqs)
    freqs <- gsub("G", " 0,0,1,0;", freqs)
    freqs <- gsub("T", " 0,0,0,1;", freqs)
    freqs <- gsub("(^[[:space:]]+|[[:space:]]+$)", "", freqs)
    freqs <- substr(freqs, 1, nchar(freqs)-1)
    root_freqs <- paste0('<rootfreqseq id="seq_Root', id, "_", clone@clone,
     '" spec="Sequence" taxon="Root', 
      id, "_", clone@clone, '" uncertain="true"
    totalcount="4" value="', freqs, '"/>', sep="")
    return(root_freqs)
  }
  root_freqs <- paste0('<rootfreqseq id="seq_Root" spec="Sequence" taxon="Root"
      totalcount="4" value="', clone@germline,'"/>', sep="")
  return(root_freqs)
}

#' Takes an airr clone object and returns BEAST2 XML for MRCA prior of the observed sequences
#' 
#' @param    clone                    an \code{airrClone} object
#' @param    id                       unique identifer for this analysis
#'
#' @return   String of XML setting the MRCA prior of the observed sequences
#'  
create_MRCA_prior_observed <- function(clone, id) {
  taxa <- paste0('<taxon id="', clone@data$sequence_id, '" spec="Taxon"/>', collapse="\n")
  distribution_xml <- 
    paste0('<distribution id="obs.prior" spec="beast.base.evolution.tree.MRCAPrior" monophyletic="true" tree="@Tree.t:', 
      id, "_", clone@clone, '">\n', 
           '<taxonset id="obs" spec="TaxonSet">\n', 
           taxa, 
           '\n</taxonset>\n',
           '</distribution>', sep="")
  return(distribution_xml)
}

#' Takes an airr clone object and returns BEAST2 XML for MRCA prior of the germline sequence
#' 
#' @param    clone                    an \code{airrClone} object
#' @param    id                       unique identifer for this analysis
#' @param    germline_range           Possible date range of germline tip
#'
#' @return   String of XML setting the MRCA prior of the germline sequence
#'  
create_MRCA_prior_germline <- function(clone, id, germline_range) {
  if(length(germline_range) != 2){
    stop("germline_range must be a vector of length 2")
  }
  taxa <- paste0('<taxon id="', 'Germline', '" spec="Taxon"/>', collapse="\n")
  distribution_xml <- 
    paste0('<distribution id="germ1.prior" spec="beast.base.evolution.tree.MRCAPrior" tipsonly="true" tree="@Tree.t:',
     id, "_", clone@clone, '">\n', '<taxonset id="germSet" spec="TaxonSet">\n', 
           taxa, 
           '\n</taxonset>\n',
           '<Uniform id="Uniform.1:germ" name="distr" lower = "',germline_range[1],'" upper="',germline_range[2],'"/>\n',
           '</distribution>', sep="")
  return(distribution_xml)
}

#' Takes an airr clone object and returns BEAST2 XML to set a height prior
#' 
#' @param    clone                    an \code{airrClone} object
#' @param    id                       unique identifer for this analysis
#' @param    start_date               starting date to use as prior, in forward time
#'
#' @return   String of XML setting the height prior
#'  
create_height_prior <- function(clone, id, start_date) {
  # start_date is in forward time, the date not the height
  distribution_xml <- 
    paste0('<distribution id="height.prior" spec="beast.base.evolution.tree.MRCAPrior" monophyletic="true" tree="@Tree.t:', 
      id, "_", clone@clone, '">\n', 
           paste0('<taxonset idref="TaxonSet.',paste0(id, "_", clone@clone),'" spec="TaxonSet">\n'), 
           '</taxonset>\n',
           '<LaplaceDistribution id="LaplaceDistribution.1" name="distr">\n',
            paste0('<parameter id="RealParameter.32" spec="parameter.RealParameter" estimate="false" name="mu">',start_date,'</parameter>\n'),
            '<parameter id="RealParameter.33" spec="parameter.RealParameter" estimate="false" name="scale">0.001</parameter>\n',
            '</LaplaceDistribution>\n',
           '</distribution>', sep="")
  return(distribution_xml)
}

#' Takes an airr clone object and returns BEAST2 XML to set a maximum height prior
#' 
#' @param    clone                    an \code{airrClone} object
#' @param    id                       unique identifer for this analysis
#' @param    max_start_date           max start date to use for prior, in forward time
#'
#' @return   String of XML setting the MRCA prior of the observed sequences
#'  
create_max_height_prior <- function(clone, id, max_start_date) {
  # max_start_date is in forward time, the date not the height
  distribution_xml <- 
    paste0('<distribution id="height.prior" spec="beast.base.evolution.tree.MRCAPrior" monophyletic="true" tree="@Tree.t:', 
      id, "_", clone@clone, '">\n', 
           paste0('<taxonset idref="TaxonSet.',paste0(id, "_", clone@clone),'" spec="TaxonSet">\n'), 
           '</taxonset>\n',
           '<Uniform id="Uniform.root" name="distr" lower="',max_start_date,'" upper="Infinity"/>',
            '</distribution>', sep="")
  return(distribution_xml)
}

#' Takes an airr clone object and returns BEAST2 XML for a trait/traitSet from a column
#' 
#' @param    clone                    an \code{airrClone} object
#' @param    trait_name               name of the trait
#' @param    column                   column in the clone data to use for the trait
#' @param    id                       unique identifer for this analysis
#' @param    trait_data_type          optional data type for the trait
#' @param    isSet                    is this a traitSet (TRUE) or a trait (FALSE)?
#' @param    include_germline_as_tip  include the germline as a tip
#'
#' @return   String of XML of the trait or traitSet
#'  
create_traitset <- function(clone, trait_name, column, id, trait_data_type=NULL, 
  isSet=FALSE, include_germline_as_tip=FALSE) {

  all_traits <- paste(clone@data$sequence_id, clone@data[[column]], 
    collapse=",\n", sep="=")
  if (include_germline_as_tip) {
    all_traits <- paste(all_traits, paste0('Germline','=', '?'), sep=",\n")
  }
  tagname <- "trait" 
  if (isSet) {
    tagname <- "traitSet"
  }
  traitset_xml <- 
    paste0('<', tagname ,' id="', trait_name, ":", id, "_", clone@clone, 
           '" spec="beast.base.evolution.tree.TraitSet" taxa="@TaxonSet.', 
           id, "_", clone@clone, 
           '" traitname="', trait_name,
           '">\n', 
           all_traits, 
           '</',tagname,'>\n',
           trait_data_type)
  
  return(traitset_xml)
}

#' Takes an airr clone object and tree and returns BEAST2 XML for setting the starting tree
#' 
#' @param    clone                    an \code{airrClone} object
#' @param    id                       unique identifer for this analysis
#' @param    tree                     starting tree, either a phylo object or a newick string
#' @param    include_germline_as_tip  include the germline as a tip
#' @param    tree_states              use states in the starting tree?
#' @param    start_edge_length        edge length to use for all branches in starting tree 
#'
#' @return   String of XML setting the starting tree
#'  
create_starting_tree <- function(clone, id, tree, include_germline_as_tip, tree_states, start_edge_length) {
  # create a starting tree in newick format
  if (inherits(tree, "phylo")) {
    ntips = length(tree$tip.label)
    if(tree_states){
      tree$node.label = (ntips + 1): (ntips + 1 + tree$Nnode - 1)
      tree$edge.length <- rep(start_edge_length, length(tree$edge.length)) # set all edge lengths to 10
    }else{
      tree$node.label = NULL
      #tree$edge.length = rep(start_edge_length, length(tree$edge.length))
      tree <- ape::multi2di(tree)
      terminal <- tree$edge[,2]<= length(tree$tip.label)
      tree$edge.length[terminal] <- start_edge_length
      tree$edge.length[!terminal] <- 1
    }
    if (!include_germline_as_tip) {
      # remove the germline tip if it exists but after numbering the nodes
      # so that the node numbers are correct
      tree <- ape::drop.tip(tree, "Germline")
    }
    newick <- ape::write.tree(tree)
  } else if (inherits(tree, "character")) {
    newick <- tree
  } else {
    stop("tree must be a phylo object or a character string in newick format")
  }
  
  # create the starting tree XML
  if(tree_states){
    starting_tree_xml <- 
    paste0('<init spec="beast.base.evolution.tree.TreeParser" id="NewickTree.t:', 
           id, "_", clone@clone, '" initial="@Tree.t:', 
           id, "_", clone@clone, '" taxa="@', 
           id, "_", clone@clone, '" IsLabelledNewick="false" adjustTipHeights="true" newick="', 
           newick, '"/>')
  }else{
    starting_tree_xml <- 
      paste0('<init spec="beast.base.evolution.tree.TreeParser" id="NewickTree.t:', 
           id, "_", clone@clone, '" initial="@Tree.t:', 
           id, "_", clone@clone, '" taxa="@', 
           id, "_", clone@clone, '" IsLabelledNewick="false" adjustTipHeights="true" newick="', 
           newick, '"/>')
  }
  
  return(starting_tree_xml)
}


#' Takes an airr clone object and template and writes a BEAST2 XML file 
#' 
#' @param    clone                    an \code{airrClone} object
#' @param    file                     output file path
#' @param    id                       unique identifer for this analysis
#' @param    time                     name of column representing sample time
#' @param    trait                    name of column representing a trait
#' @param    trait_data_type          optional data type for the trait
#' @param    template                 XML template
#' @param    mcmc_length              number of MCMC iterations
#' @param    log_every                frequency of states logged. \code{auto} will divide mcmc_length by log_target
#' @param    replacements             list of additional replacements to make in the template
#' @param    include_germline_as_root include germline in analysis as root?     
#' @param    include_germline_as_tip  include germline in analysis as tip?     
#' @param    germline_range           possible date range of germline
#' @param    tree                     starting tree, either a phylo object or a newick string
#' @param    trait_list               list of all possible trait values
#' @param    log_every_trait          frequency of trait states logged relative to log_every
#' @param    tree_states              use states in the starting tree?
#' @param    start_edge_length        edge length to use for all branches in starting tree
#' @param    start_date               starting date to use as prior, in forward time
#' @param    max_start_date           max starting date to use as prior, in forward time
#' @param    ...                      additional arguments for XML writing functions
#'
#' @return   File path of the written XML file
#'  
write_clone_to_xml <- function(clone, file, id, time=NULL, trait=NULL, 
  trait_data_type=NULL, template=NULL, mcmc_length=1000000, log_every=1000, replacements=NULL, 
  include_germline_as_root=FALSE, include_germline_as_tip=FALSE, 
  germline_range=c(-10000,10000), tree=NULL, trait_list=NULL, log_every_trait=10, tree_states=FALSE,
  start_edge_length=100, start_date=NULL, max_start_date=NULL,...) {
  
  kwargs <- list(...)

  # read in a template file
  if (is.null(template)) {
    # template <- system.file("extdata", "template.xml", package = "scoper")
    template = "template.xml"
  }
  xml <- readLines(template)

  xml <- gsub('chainLength="[^"]*"', 
    paste0('chainLength="', format(mcmc_length, scientific=F), '"', sep=''), xml) 
  
  # log traited tree 10x less frequently than regular tree
  trait_logger_index <- grepl("treeWithTraitLogger", xml)
  xml[trait_logger_index] <- gsub('logEvery="[^"]*"', 
    paste0('logEvery="', format(log_every*log_every_trait, scientific=F), '"', sep=''), 
    xml[trait_logger_index])
  
  xml[!trait_logger_index] <- gsub('logEvery="[^"]*"', 
    paste0('logEvery="', format(log_every, scientific=F), '"', sep=''), 
    xml[!trait_logger_index])
  
  data <- create_alignment(clone, id, include_germline_as_tip)
  # replace the ${DATA} placeholder with actual data
  xml <- gsub("\\$\\{DATA\\}", data, xml)
  xml <- gsub("\\$\\{CLONE\\}", paste0(id, "_", clone@clone), xml)
  
  # if the date argument is null but ${DATE} is in the template file, raise an error
  if (is.null(time) && any(grepl("\\$\\{DATE\\}", xml))) {
    stop("Date argument is NULL but ${DATE} is in the template file")
  }
  if (!is.null(time)) {
    date_trait <- create_traitset(clone, "date", time, id)
    # replace the ${DATE} placeholder with the dates
    xml <- gsub("\\$\\{DATE\\}", date_trait, xml)
  }
  
  if (is.null(trait) && any(grepl("\\$\\{TRAIT\\}", xml))) {
    stop("Trait argument is NULL but ${TRAIT} is in the template file")
  }
  if (!is.null(trait) && is.null(trait_data_type)) {
    stop("Trait argument is given but no trait data type was provided")
  }
  if (!is.null(trait)) {
    sample_trait <- create_traitset(clone, "newTrait", trait, id, 
      trait_data_type, isSet=TRUE, include_germline_as_tip=include_germline_as_tip)
    # replace the ${TRAIT} placeholder with the sample trait
    xml <- gsub("\\$\\{TRAIT\\}", sample_trait, xml)
    if (any(grepl("\\$\\{TRAIT_NAME\\}", xml))) {
      xml <- gsub("\\$\\{TRAIT_NAME\\}", trait, xml)
    } else {
      # if there is no ${TRAIT_NAME} placeholder, assume these are old xml templates
      # and replace tag="location" with tag="trait"
      xml <- gsub('tag="location"', paste0('tag="',trait,'"'), xml)
    }
  }
  if (any(grepl("\\$\\{NODES\\}", xml))) {
    # replace the ${NODES} placeholder with the number of nodes in this tree
    tips <- nrow(clone@data) + 1 # node numbering includes germline as a tip
    nodes <- 2*tips-1
    xml <- gsub("\\$\\{NODES\\}", nodes, xml)
  }
  
  if (any(grepl("\\$\\{MRCA\\}", xml))) {
    # replace the ${MRCA} placeholder with the mrca prior
    mrca_priors <- ""
    if (include_germline_as_tip) {
      mrca_priors <- paste0(create_MRCA_prior_observed(clone, id), 
        create_MRCA_prior_germline(clone, id, germline_range), sep="\n")
    }
    if(!is.null(start_date)){
      mrca_priors <- paste0(mrca_priors, "\n",
        create_height_prior(clone, id, start_date) )
    }
    if(!is.null(max_start_date)){
      mrca_priors <- paste0(mrca_priors, "\n",
        create_max_height_prior(clone, id, max_start_date) )
    }
    
    xml <- gsub("\\$\\{MRCA\\}", mrca_priors, xml)
  }

  if (any(grepl("\\$\\{OPERATORS\\}", xml))) {
    # replace the ${OPERATORS} placeholder with the operators we want to add
    # can add other potential operators here
    operators <- ""
    if (include_germline_as_tip) {
      
      operators <-  
      paste0('<operator id="TipDatesRandomWalker.01" windowSize="1" spec="beast.base.evolution.operator.TipDatesRandomWalker" taxonset="@germSet" tree="@Tree.t:',
        id, '_', clone@clone,'" weight="1.0"/>\n')
    }
    
    xml <- gsub("\\$\\{OPERATORS\\}", operators, xml)
  }
  if (any(grepl("\\$\\{ROOTFREQS\\}", xml))) { 
    root_freqs <- ""
    if (include_germline_as_root) {
      # replace spec="TreeLikelihood" with spec="rootfreqs.TreeLikelihood"
      xml <- gsub('spec="TreeLikelihood"', 'spec="rootfreqs.TreeLikelihood"', xml)
      # replace the ${ROOTFREQS} placeholder with the root frequencies
      root_freqs <- create_root_freqs(clone, id)
    }
    xml <- gsub("\\$\\{ROOTFREQS\\}", root_freqs, xml)
  }
  if (any(grepl("\\$\\{EMP_EQFREQS\\}", xml))) { 
    getFreqs = function(seqs){
      nts <- c("A", "C","G","T")
      splits <- strsplit(seqs, split="")
      counts <- table(unlist(splits))[nts]
      freqs <- counts/sum(counts)
      freqs
    }
    freqs <- getFreqs(clone@data$sequence)
    freqstring <- paste(freqs, collapse=" ")
    print(freqstring)
    xml <- gsub("\\$\\{EMP_EQFREQS\\}", freqstring, xml)
  }

    if (!is.null(tree)) {
    # replace the random tree with the provided tree in newick format
    start_idx <- grep("<init", xml)
    end_idx <- grep("</init>", xml)
    if (length(start_idx) == 1 && length(end_idx) == 1 && start_idx < end_idx) {
      xml <- c(
        xml[1:(start_idx-1)],
        c(create_starting_tree(clone, id, tree, include_germline_as_tip, tree_states, start_edge_length)),
        xml[(end_idx+1):length(xml)]
      )
      # TODO: check if there are traits associated with the internal nodes
      # if so, add the integer corresponding to each node's trait to an array
      # such that starting_traits[i] is the trait for node i
      if(tree_states){
        tips <- nrow(clone@data) + 1 # node numbering includes germline as a tip
        nodes <- 2*tips-1
        if("state" %in% names(tree)){
          states <- strsplit(tree$state, split=",")
          starting_traits_all <- sapply(states, function(x)x[1])
          starting_traits <- match(starting_traits_all, trait_list)-1
          if(sum(is.na(starting_traits)) > 0){
            print(trait_list)
            stop(paste0("unknown state found", paste0(starting_traits_all, collapse=" ")))
          }
        }else{
          warning("States not found in starting tree, setting to 0")
          starting_traits <- rep(0, nodes)
        }
        # replace the value of the traitCategories parameter with the starting traits
        trait_categories_index <- grep("id='traitCategories'", xml)
        if (length(trait_categories_index) > 0) {
          xml[trait_categories_index] <- 
            gsub('value="[^"]*"', 
                 paste0('value="', paste(starting_traits, collapse=" "), '"'), 
                 xml[trait_categories_index])
        } else {
          stop("Could not find traitCategories parameter in the template file")
        }
      }
    } else {
      stop("Could not find <init> tag in the template file")
    }
  }

  matches <- unlist(regmatches(xml, gregexpr("\\$\\{([^}]+)\\}", xml)))
  template_variables <- unique(sub("\\$\\{([^}]+)\\}", "\\1", matches))
  # print(paste("Template variables found:", paste(template_variables, collapse=", ")))

  for (var in template_variables) {
    if (!var %in% names(kwargs)) {
      stop(paste("Variable", var, "not provided but is required by the template."))
    }
    if (var %in% replacements) {
      next # skip replacements, they are handled later
    }

    value <- kwargs[[var]]

    # replace the ${VAR} placeholder with the value from kwargs
    xml <- gsub(paste0("\\$\\{", var, "\\}"), value, xml)
  }

  # this is useful for e.g. GIBLE templates where we want to replace some variables later than others
  if (!is.null(replacements)) {
    for (replacement in replacements) {
      value <- kwargs[[replacement]]

      # replace the ${VAR} placeholder with the value from kwargs
      xml <- gsub(paste0("\\$\\{", replacement, "\\}"), value, xml)
    }
  }
  # open a connection to the file
  file <- paste0(file, "_", clone@clone, ".xml")
  con <- file(file, "w")
  
  # write the XML file
  writeLines(xml, con)
  
  # close the connection
  close(con)

  return(file)
}

#' Wrapper to write multiple clones to XML files 
#' 
#' @param    data                     a list of \code{airrClone} objects
#' @param    id                       identifer for this analysis
#' @param    trees                    optional list of starting trees, either phylo objects or newick strings
#' @param    time                     name of column representing sample time
#' @param    trait                    name of column representing a trait
#' @param    template                 XML template
#' @param    outfile                  output file path prefix
#' @param    replacements             list of additional replacements to make in the template
#' @param    trait_list               list of all possible trait values
#' @param    mcmc_length              number of MCMC iterations
#' @param    log_every                frequency of states logged. \code{auto} will divide mcmc_length by log_target         
#' @param    include_germline_as_root include germline in analysis as root?     
#' @param    include_germline_as_tip  include germline in analysis as tip?     
#' @param    germline_range           possible date range of germline
#' @param    tree_states              use states in the starting tree?
#' @param    start_edge_length        edge length to use for all branches in starting tree
#' @param    start_date               starting date to use as prior, in forward time
#' @param    max_start_date           max starting date to use as prior, in forward time
#' @param    ...                      additional arguments for XML writing functions
#'
#' @return   File paths of the written XML files
#'  
write_clones_to_xmls <- function(data, id, trees=NULL, time=NULL, trait=NULL, template=NULL, 
  outfile=NULL, replacements=NULL, trait_list=NULL, 
  mcmc_length=1000000, log_every=1000, include_germline_as_root=FALSE, 
  include_germline_as_tip=FALSE, germline_range=c(-10000,10000), 
  tree_states=FALSE, start_edge_length=100, start_date=NULL, max_start_date=NULL,
  ...) {

  kwargs <- list(...)

  # iterate over the clones to first create trait data type if trait exists
  if (!is.null(trait)) {
    if (is.null(trait_list)) {
      # get all the possible values of the trait
      traits <- c()
      for (i in 1:length(data)) {
        traits <- c(traits, unique(data[[i]]@data[[trait]]))
      }
      trait_list <- sort(unique(traits))
    }
    codeMap <- paste(trait_list, 0:(length(trait_list)-1), collapse=",\n", sep="=")
    codeMap <- paste0(codeMap, ",\n? = ", paste0(0:(length(trait_list)-1), collapse=" "))
    trait_data_type <- paste0('<userDataType id="traitDataType.newTrait" spec="beast.base.evolution.datatype.UserDataType" codeMap="',
     codeMap, '" codelength="-1" states="', length(trait_list), '"/>')
  }
  xmls = c()
  for (i in 1:length(data)){
    #if ("trees" %in% names(kwargs)) {
    if(!is.null(trees)){
      # if trees are provided, use them
      tree <- trees[[i]]
    } else {
      tree <- NULL
    }
    xmls = c(xmls, write_clone_to_xml(data[[i]], 
                     file=outfile, 
                     id=id, 
                     time=time, 
                     trait=trait, 
                     trait_data_type=trait_data_type, 
                     template=template,
                     replacements=replacements, 
                     include_germline_as_root=include_germline_as_root,
                     include_germline_as_tip=include_germline_as_tip,
                     mcmc_length=mcmc_length,
                     log_every=log_every,
                     germline_range=germline_range,
                     tree=tree,
                     trait_list=trait_list,
                     tree_states=tree_states,
                     start_date=start_date,
                     max_start_date=max_start_date,
                     ...))
  }
  return(xmls)
}


#' Reads in a BEAST output directory
#' 
#' \code{readBEAST} Reads in data from BEAST output directory
#' @param clones     either a tibble (getTrees) or list of \code{airrClone} object
#' @param beast      location of beast binary directory (beast/bin)
#' @param dir        directory where temporary files will be placed.
#' @param id         unique identifer for this analysis
#' @param trait      Trait coolumn used         
#' @param asr        Log ancestral sequences?
#' @param full_posterior  Read un full distribution of parameters and trees?
#' @param nproc      Number of cores for parallelization. Uses 1 core per tree.
#' @param quiet      amount of rubbish to print to console
#' @param burnin         percent of initial tree samples to discard (1-100)
#' @param nproc          cores to use
#' @param low_ram        run with less memory (slower)       
#'
#' @return   
#' If data is a tibble, then the input clones tibble with additional columns for 
#' trees and parameter estimates given the specified burnin. If input is just a 
#' list of airrClone objects, it will return the corresponding list of trees
#' given the burnin
#'  
#' @export
readBEAST <- function(clones, dir, id, beast, burnin=10, trait=NULL, nproc = 1, 
  quiet=0, full_posterior=FALSE, asr=FALSE, low_ram=TRUE) {

  if(!"list" %in% class(clones) && "data" %in% names(clones)){
    data <- clones$data
  }else if("list" %in% class(clones)){
    data <- clones
  }else{
    stop("Input data type not supported")
  }

  beast <- path.expand(beast)
  
  annotator_exec <- file.path(beast,"treeannotator")
  if(file.access(annotator_exec, mode=1) == -1) {
    stop("The file ", annotator_exec, " cannot be executed.")
  }
  analyser_exec <- file.path(beast,"loganalyser")
  if(file.access(annotator_exec, mode=1) == -1) {
    stop("The file ", annotator_exec, " cannot be executed.")
  }
 # Run treeannotator in parallel
  capture <- parallel::mclapply(1:length(data), function(x) {
    y <- data[[x]]@clone
    if(is.null(trait)){
      treesfile <- ifelse(asr, paste0(id, "_", data[[x]]@clone, "_asr.trees"),
                      paste0(id, "_", data[[x]]@clone, ".trees"))
      treesfile <- file.path(dir, treesfile)
      treefile <- file.path(dir, paste0(id, "_", data[[x]]@clone, ".tree"))
    }else{
      treesfile <- ifelse(asr, paste0(id, "_", data[[x]]@clone, "_asr.trees"),
                      paste0(id, "_", data[[x]]@clone, "_tree_with_trait.trees"))
      treesfile <- file.path(dir, treesfile)
      treefile <- file.path(dir, paste0(id, "_", data[[x]]@clone, 
        "_tree_with_trait.tree"))
    }
    command <- paste("-burnin", burnin, treesfile, treefile)
    if(low_ram){
      command <- paste("-lowMem TRUE", command)
    }

    console_log <- file.path(dir, paste0(id, "_", data[[x]]@clone,".log"))

    #command <- paste(command, ">>", console_log)

    if(quiet < 1){
      print(paste(annotator_exec,command))
    }
      
    params <- list(annotator_exec, command, stdout=TRUE, stderr=TRUE)
  
    status <- tryCatch(do.call(base::system2, params), error=function(e){
         print(paste("TreeAnnotator error: ",e));
         return(e)
     }, warning=function(w){
         print(paste("TreeAnnotator warnings ",w));
         return(w)
     })

    status
    }, mc.cores=nproc)

  for(i in 1:length(capture)){
      if("error" %in% class(capture[[i]]) || grepl("[E|e]rror",capture[i])){
        print(capture[[i]])
        stop(paste("Error running Treeannotator (see above), clone", data[[i]]@clone))
      }
    }

  # Run loganalyser in parallel
  capture <- parallel::mclapply(1:length(data), function(x) {
    y <- data[[x]]@clone
    logfile <- file.path(dir, paste0(id, "_", data[[x]]@clone,".log"))
    outfile <- file.path(dir, paste0(id, "_", data[[x]]@clone,"_log.tsv"))
    command <- paste("-quiet -b", burnin, logfile, ">", outfile)
    
    if(quiet < 1){
      print(paste(analyser_exec,command))
    }
      
    params <- list(analyser_exec, command, stdout=TRUE, stderr=TRUE)
  
    status <- tryCatch(do.call(base::system2, params), error=function(e){
         print(paste("Loganalyser error: ",e));
         return(e)
     }, warning=function(w){
         print(paste("Loganalyser warnings ",w));
         return(w)
     })

    status
    }, mc.cores=nproc)


  for(i in 1:length(capture)){
    if("error" %in% class(capture[[i]]) || grepl("[E|e]rror",capture[i])){
      print(capture[[i]])
      stop(paste("Error running Loganalyser (see above), clone", data[[i]]@clone))
    }
  }
  
  # read in tree and parameter log
  # TODO add nodes and sequences to tree
  trees <- list()
  for(i in 1:length(data)){
    if(is.null(trait)){
      treefile <- file.path(dir, paste0(id,"_",data[[i]]@clone, ".tree"))
    }else{
      treefile <- file.path(dir, paste0(id,"_",data[[i]]@clone, "_tree_with_trait.tree"))
    }
    logfile <- file.path(dir, paste0(id,"_",data[[i]]@clone, ".log"))
    logoutfile <- file.path(dir, paste0(id, "_", data[[i]]@clone,"_log.tsv"))

    beast <- treeio::read.beast(treefile)
    if("error" %in% class(beast)){
      stop(paste("Couldn't read in ",treefile))
    }
    l <- readLines(logoutfile)
    log <- read.table(text=l[4:(length(l)-1)], header=TRUE)

    # add parameter summary
    beast@info$parameters <- log
    if(full_posterior){ 
      treesfile <- file.path(dir, paste0(id,"_",data[[i]]@clone, ".trees"))
      l <- readLines(treesfile, warn=FALSE)
      if(!grepl("End;",l[length(l)])){
        l[length(l) + 1] = "End;"
        warning("Adding End; to ",treesfile)
        #make new file to avoid overwriting
        treesfile <- file.path(dir, paste0(data[[i]]@clone, "_end.trees"))
        writeLines(l, con=treesfile)
      }
      phylos <- ape::read.nexus(treesfile)
      burn <- floor(length(phylos)*burnin/100)
      phylos <- phylos[(burn+1):length(phylos)]
      phylos <- lapply(phylos, function(y){
        y$tip.label <- sapply(strsplit(y$tip.label,"_"), function(x)
          paste0(x[1:(length(x)-1)], collapse="_"))
        y
      })
      
      beast@info$tree_posterior <- phylos

      l <- read.table(logfile, header=TRUE)
      beast@info$parameters_posterior <- tidyr::gather(l, "parameter", "value", -(!!rlang::sym("Sample")))
    }
    beast@info$name <- data[[i]]@clone
    trees[[i]] <- beast
  }

  if(quiet < 1)print("Ran readBEAST")

  if(!"list" %in% class(clones) && "data" %in% names(clones)){
    clones$trees <- trees
    clones$parameters <- lapply(trees, function(x)x@info$parameters)
    #clones <- getParams(clones, burnin=burnin)
    return(clones)
  }else if("list" %in% class(clones)){
    return(trees)
  }else{
    stop("Input data type not supported")
  }
}

#' get values for Bayesian Skyline plot
#' 
#' \code{makeSkyline} 
#' @param  logfile   Beast log file
#' @param  treesfile BEAST trees file 
#' @param  burnin    Burnin percentage (1-100) 
#' @param  bins      number of bins for plotting
#' @param  youngest  timepoint of the most recently tip sampled (if 0, backward time used)
#' @param  clone_id  name of the clone being analyzed (if desired)
#' @param  max_height max height to use (min, median, mean, max)
#' @return   Bayesian Skyline values for given clone
#'
#' @export
makeSkyline <- function(logfile, treesfile, burnin, bins=100, youngest=0, 
    clone_id=NULL, max_height=c("min","median","mean","max")){
    
    l <- tryCatch(read.csv(logfile, header=TRUE, sep="\t", comment.char="#"),error=function(e)e)
    if("error" %in% class(l)){
        stop(paste("couldn't open",logfile))
    }
    phylos <- tryCatch(ape::read.nexus(treesfile), error=function(e)e)
    if("error" %in% class(phylos)){
        stop(paste("couldn't open", treesfile))
    }
    params <- tidyr::gather(l, "parameter", "value", -(!!rlang::sym("Sample")))

    if(!"bPopSizes.1" %in% unique(params$parameter)){
        stop(paste("log file doesn't have pop sizes.",
            "Was it run with skyline tree_prior='coalescent_skyline'?"))
    }

    burn <- floor(length(phylos)*burnin/100)
    samples <- unique(params$Sample)
    if(burn > 0){
      phylos <- phylos[(burn+1):length(phylos)]
      samples <- samples[(burn+1):length(samples)]
      params <- dplyr::filter(params, !!rlang::sym("Sample") %in% samples)
    }
    if(dplyr::n_distinct(params$Sample) != length(phylos)){
      warning("Parameter and tree posteriors not same length, subsetting")
      treestates <- as.numeric(gsub("STATE_","",names(phylos)))
      commonstates <- intersect(treestates, params$Sample)
      phylos <- phylos[treestates %in% commonstates]
      params <- dplyr::filter(params, !!rlang::sym("Sample") %in% commonstates)
      samples <- unique(params$Sample)
    }

    groups <- dplyr::filter(params, grepl("GroupSizes", !!rlang::sym("parameter")))
    pops <- dplyr::filter(params, grepl("PopSizes", !!rlang::sym("parameter")))

    if(sum(pops$value < 0) > 0){
        stop(paste(logfile, "found popsizes < 0, can't continue"))
    }
    if(sum(groups$value < 0) > 0){
        stop(paste(logfile, "found groupsizes < 0, can't continue"))
    }

    pops$index <- as.numeric(gsub("bPopSizes\\.","",pops$parameter))
    groups$index <- as.numeric(gsub("bGroupSizes\\.","",groups$parameter))

    # smallest tree height in log file (this is what tracer seems to do)
    if(max_height == "min"){
      maxheight <- min(dplyr::filter(params, !!rlang::sym("parameter") == "TreeHeight")$value)
    }else if(max_height == "median"){
      maxheight <- stats::median(dplyr::filter(params, !!rlang::sym("parameter") == "TreeHeight")$value)
    }else if(max_height == "mean"){
      maxheight <- mean(dplyr::filter(params, !!rlang::sym("parameter") == "TreeHeight")$value)
    }else if(max_height == "max"){
      maxheight <- max(dplyr::filter(params, !!rlang::sym("parameter") == "TreeHeight")$value)
    }else{
      stop("max_height option must be min, median, mean, or max")
    }
    
    if(youngest > 0){
        mintime <- youngest - maxheight
        maxtime <- youngest
    }else{
        mintime <- 0
        maxtime <- maxheight - youngest
    }

    binwidth <- (maxtime - mintime)/(bins - 1)

    all_intervals <- dplyr::tibble()
    for(index in 1:length(phylos)){
      tr <- phylos[[index]]
      sample <- samples[index]
      mrca <- ape::getMRCA(tr, tip=tr$tip.label)
      d <- ape::dist.nodes(tr)
      times <- d[mrca,]
      maxheight <- max(times)
      nodes <- maxheight - times[(length(tr$tip.label)+1):length(times)]
      nodes <- nodes[order(nodes, decreasing=FALSE)]


      pop <- dplyr::filter(pops, !!rlang::sym("Sample") == sample)
      group <- dplyr::filter(groups, !!rlang::sym("Sample") == sample)

      temp <- nodes
      results <- dplyr::tibble()
      for(i in 1:nrow(group)){
        groupsize <- group$value[i]
        popsize <- pop$value[i]
        events <- temp[1:(groupsize)]
        results <- bind_rows(results,
          dplyr::tibble(end=events[length(events)], interval=i,
            events=length(events), popsize=popsize))
        temp <- temp[-1:-(groupsize)]
      }
      results$sample <- sample
      results$index <- index
      all_intervals <- bind_rows(all_intervals, results)
    }

    indistinct <- all_intervals %>%
        group_by(sample) %>%
        summarize(distinct = dplyr::n_distinct(!!rlang::sym("end")),
            n = n()) %>%
        dplyr::filter(!!rlang::sym("distinct") < n) %>%
        pull(sample)

    if(length(indistinct) > 0){
        warning(paste(logfile, "Removing",length(indistinct),
            "samples with indistinct intervals. This shouldn't happen."))
        all_intervals <- dplyr::filter(all_intervals, !(!!rlang::sym("sample") %in% indistinct))
        if(nrow(all_intervals) == 0){
            stop("No intervals left :-(")
        }
    }
 
    skyline <- tidyr::tibble()
    n_sample <- dplyr::n_distinct(all_intervals$sample)
    bin_intervals <- seq(0, length=bins, by=binwidth)
    interval_bins <- dplyr::tibble()
    for(j in 1:(length(bin_intervals)-1)){
      bin_start <- bin_intervals[j]
      bin_end <- bin_intervals[j+1]

      matches <- all_intervals %>%
        dplyr::group_by(!!rlang::sym("sample")) %>%
        dplyr::filter(!!rlang::sym("end") > bin_start) %>%
        dplyr::slice_min(!!rlang::sym("end"))

#      if(nrow(matches) != n_sample){
#        stop("didn't find some indexes")
#      }

      matches <- select(matches, !!rlang::sym("end"), 
        !!rlang::sym("popsize"), !!rlang::sym("sample"), !!rlang::sym("index"))
      matches$bin <- bin_start
      skyline <- bind_rows(skyline, matches)
    }

    # similar behavior to 
    # dr.stat.getQuantile in beast-mcmc
    getQuantiles <- function(x, probs=NULL){
      if(is.null(probs)){
        stop("probs can't be null")
      }
      x <- sort(x)
      index <- sapply(probs, function(y)
        ceiling(length(x)*y))
      x <- x[index]
      names(x) <- paste0(probs*100,"%")
      return(x)
    }

    skyplot <- skyline %>%
      dplyr::group_by(!!rlang::sym("bin")) %>%
      dplyr::summarize(
        mean=mean(!!rlang::sym("popsize")), 
        median=getQuantiles(!!rlang::sym("popsize"), 0.5),
        lci=getQuantiles(!!rlang::sym("popsize"), 0.025), 
        uci=getQuantiles(!!rlang::sym("popsize"), 0.975))

    if(!is.null(clone_id)){
        skyplot$clone_id <- clone_id
    }
    if(youngest != 0){
        skyplot$bin <- youngest - skyplot$bin
    }

    return(skyplot)
}

#' Make data frames for Bayesian skyline plots
#' 
#' \code{makeSkylines} 
#' @param  clones    clone tibble
#' @param  dir       directory of BEAST trees file 
#' @param  id        unique identifer for this analysis
#' @param  time      name of time column
#' @param  bins      number of bins for plotting
#' @param  burnin    Burnin percent (default 10) 
#' @param  verbose   if 1, print name of clones
#' @param  forward   plot in forward or (FALSE) backward time?
#' @param  nproc     processors for parallelization (by clone)
#' @param  max_height max height to use (min, median, mean, max)
#' @return   Bayesian Skyline values for given clone
#' @details Burnin set from readBEAST or getTrees
#' @export
getSkylines <- function(clones, dir, id, time, burnin=10, bins=100, verbose=0, forward=TRUE,
    nproc=1, max_height=c("min","median","mean","max")){

    treesfiles <- sapply(clones$data, function(x)
        file.path(dir, paste0(id, "_", x@clone, ".trees")))

    logfiles <- sapply(clones$data, function(x)
        file.path(dir, paste0(id, "_", x@clone, ".log")))

    if(forward){
        youngest <- sapply(clones$data, function(x)
            max(as.numeric(x@data[[time]])))
    }else{
        youngest <- rep(0, length=nrow(clones))
    }

    skylines <- parallel::mclapply(1:nrow(clones), function(x){
        if(verbose != 0){
            print(paste(clones$clone_id[x], logfiles[x], 
                treesfiles[x], youngest[x]))
        }
        tryCatch(makeSkyline(logfile=logfiles[x], treesfile=treesfiles[x],
            youngest=youngest[x], burnin=burnin, bins=bins, 
            clone_id=clones$clone_id[x], max_height=max_height), error=function(e)e)
    }, mc.cores=nproc)

    clones$skyline <- skylines
    for(i in 1:nrow(clones)){
        if("error" %in% class(skylines[[i]])){
            print(skylines[[i]])
            warning(paste("Error making skyline clone,",clones$clone_id[i]))
            clones$skyline[[i]] <- NA
        }
    }
    return(clones)
}