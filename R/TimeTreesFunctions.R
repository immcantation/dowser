# Functions for performing time tree analysis on B cell lineage trees

#' Iteratively resume getTimeTrees until convergence, as defined by 
#' all parameters (except those in \code{ignore} vector) having ESS 
#' greater than or equal to the specified ess_cutoff
#'
#' \code{getTimeTreesIterate} Iteratively resume getTimeTrees til convergence.
#' @param    clones     a tibble of \code{airrClone} objects, the output of
#'                      \link{formatClones}
#' @param    template   XML template
#' @param    beast      location of beast binary directory (beast/bin)
#' @param    dir        directory where temporary files will be placed.
#' @param    id         unique identifer for this analysis
#' @param    time       Name of sample time column        
#' @param    burnin     Burnin percent (default 10)                 
#' @param    trait      Trait column to be used. "?" is treated as 
#'                      ambiguous (any type) and is not counted as a type.         
#' @param    nproc      Number of processes for parallelization. At most 1 process/tree can be used. Use nproc=1 for Windows.
#' @param    posterior  Read un full distribution of parameters and trees? Can be "none" to just have
#'                      summary objects, "all" to have parameters, trees, and trees_with_traits, or a vector with the desired
#'                      combination of "parameters", "trees_with_traits", and "trees".
#' @param    quiet      amount of rubbish to print to console
#' @param    iterations Maximum number of times to run MCMC chain for each clone. The first run counts as 1.
#' @param    ess_cutoff Minimum number of ESS for all parameters
#' @param    ignore     Vector of parameters to ignore for ESS calculation
#' @param    quiet      quiet notifications if > 0
#' @param    continue   If TRUE, continue an earlier run of the same \code{id} in \code{dir}
#'                      from saved checkpoints. (default FALSE) Warning: if the settings or
#'                      the data are not the same as in the earlier run, files may be overwritten, as if continue=FALSE.
#'                      See details. 
#' @param    asr          Log ancestral sequences?
#' @param    low_ram    run with less memory (slightly slower)  
#' @param    trim_ids    remove last _ group from tips?      
#' @param    skip_errors If TRUE, a clone that fails is dropped with a warning and the other clones
#'                      continue. If FALSE (default), any failing clone stops the analysis.
#' @param    backup      If TRUE, a backup folder with the existing files of this analysis (BEAST files and
#'                      checkpoint folder) is made before the run starts. The files are moved if \code{continue}
#'                      is FALSE and copied if it is TRUE. See details. (default FALSE)
#' @param    trait_list  Optional vector of all types of the trait. If NULL, it is made from the values
#'                      in all clones, so that every clone uses the same types.
#' @param    ...        Additional arguments for getTimeTrees. If an argument is a function, it is
#'                      called with the one-row tibble of each clone, and the value that it returns
#'                      is used for that clone. A list in \code{trees} is cut to the tree of each clone.
#'
#' @details
#' 
#' For a full list of options (of which there are many more), see \code{getTimeTrees}
#' 
#' Each clone runs as its own BEAST job. A clone with ESS below \code{ess_cutoff} is resumed (its MCMC
#' chain continues) until it converges or has run \code{iterations} times. For resumability in case of interruption, some intermediate and checkpoint files are saved. 
#' 
#' \strong{continue=TRUE} is for resuming a job that was cancelled, or for running more \code{iterations}.
#' Repeat the original call, with the same \code{clones}, \code{dir}, \code{id} and settings, and add
#' \code{continue=TRUE} (with a larger \code{iterations} to run more). Clones with a checkpoint resume from
#' it, unless they already converged or reached \code{iterations}. Clones without one start from the
#' beginning.
#' You can also change \code{ess_cutoff} or \code{ignore}: convergence is checked again from the last
#' checkpoint of each clone, and clones that are now below the cutoff run again if they have run fewer than
#' \code{iterations} times.
#' 
#' Do not use \code{continue=TRUE} when changing settings (e.g. \code{seed}, \code{mcmc_length} or
#' \code{template}) or the data of an analysis that you want to keep. Each checkpoint has the XML that was
#' used for the clone, and the seed, the java flag and the BEAST path. If any of them is not the same in this
#' call, the checkpoint and the BEAST files of that clone are deleted, and the clone starts again from round 1.
#' The other clones continue, and a message lists the clones. Use a new \code{id} or \code{dir} instead. Do not
#' use if you have edited the BEAST files in \code{dir}.
#' Checkpoints from older versions of dowser cannot be used.
#' 
#' \strong{backup=TRUE} makes the folder \code{backup_<id>_<date>_<time>} in \code{dir} before the run starts, with the
#' BEAST and checkpoint files of \code{id}. Use it if your job may be
#' restarted automatically (for example by a job scheduler), e.g. the restarted job starts a new run, you cancel it and
#' use \code{\link{restoreIterateBackup}} to put the files back, and then call the function again with
#' \code{continue=TRUE}. Use it also for insurance if \code{dir} may not be empty before this call.
#' 
#' If a clone fails and \code{skip_errors=FALSE}, running clones finish their current round, no new clone
#' starts, and then the function stops with the error. 
#' 
#' If you interrupt R (for example with Ctrl-C), BEAST processes that the workers started can keep running.
#' Check for them with the appropriate tools for your system.
#' 
#' @return   A tibble of \code{tidytree} and \code{airrClone} objects. Clones that failed with \code{skip_errors=TRUE} are not in the tibble.
#'
#' @details
#' For examples and vignettes, see https://dowser.readthedocs.io
#'
#' @export
getTimeTreesIterate <- function(clones, template, beast, dir, id, time,
  iterations=10, ess_cutoff=200, burnin=10,
  ignore = c("traitfrequencies"), continue=FALSE, 
  posterior=c("none","all","parameters","trees_with_traits","trees"), 
  trait=NULL, nproc=1, quiet=0, asr=FALSE, low_ram=TRUE, trim_ids=FALSE,
  skip_errors=FALSE, trait_list=NULL, backup=FALSE, ...){

  if(is.null(dir)){
    stop("dir must be specified when running BEAST")
  }
  dir <- path.expand(dir)
  if(is.null(id)){
    id <- "sample"
  }

  dots <- list(...)
  all_ids <- as.character(clones$clone_id)
  # make the backup first
  if (isTRUE(backup)) {
    makeBackup(dir, id, all_ids, continue)
  }

  # one list of types for all clones, so that every clone has the same codes
  if(!is.null(trait) && is.null(trait_list) && nrow(clones) > 0){
    trait_list <- makeTraitList(lapply(clones$data, function(x)x@data[[trait]]))
  }

  # one round of one clone. The engine saves the checkpoint afterwards.
  work <- function(cs, resume_id){
    clone_dots <- iterateCloneDots(dots, cs, all_ids)
    do.call(getTimeTrees, c(list(clones=cs, resume_clones=resume_id, 
      posterior="none", template=template, dir=dir, id=id, beast=beast, 
      burnin=burnin, time=time, trait=trait, trait_list=trait_list, 
      quiet=quiet, nproc=1, asr=asr, trim_ids=trim_ids, low_ram=low_ram), 
      clone_dots))
  }

  # record options that are not in the xml, for the checkpoint manifest
  clone_options <- function(cs){
    clone_dots <- iterateCloneDots(dots, cs, all_ids)
    value <- function(name){
      x <- if(name %in% names(clone_dots)) clone_dots[[name]] else eval(formals(getTimeTrees)[[name]])
      if(is.numeric(x)) x <- as.numeric(x)   # 1e5 and 100000L are the same
      paste(x, collapse=",")
    }
    paste0("seed=", value("seed"), " java=", value("java"),
      " beast=", normalizePath(path.expand(beast), mustWork=FALSE))
  }

  # The md5 sum of the XML that dowser would use for a clone. BEAST is not
  # started: the function xml_check stops getTimeTrees after the XML is
  # written. The XML is written in a temporary folder. NA if the XML cannot be made.
  xml_md5 <- function(cid){
    cs <- clones[match(cid, all_ids), ]
    tmp <- tempfile("iterate_xml_")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive=TRUE))
    stop_after_xml <- function(paths){
      stop(structure(class=c("iterate_xml_only", "error", "condition"),
        list(message="XML written", call=NULL, paths=paths)))
    }
    paths <- NULL
    utils::capture.output(paths <- tryCatch(
      do.call(getTimeTrees, c(list(clones=cs, resume_clones=NULL, 
        posterior="none", template=template, dir=tmp, id=id, beast=beast, 
        burnin=burnin, time=time, trait=trait, trait_list=trait_list, 
        quiet=1, nproc=1, asr=asr, trim_ids=trim_ids, low_ram=low_ram, 
        xml_check=stop_after_xml), iterateCloneDots(dots, cs, all_ids))),
      iterate_xml_only=function(cnd)cnd$paths, error=function(e)NULL))
    if(!is.character(paths) || length(paths) == 0){
      return(NA_character_)
    }
    unname(tools::md5sum(paths[1]))
  }

  # With continue=TRUE, a clone whose checkpoint has other options or another XML starts again
  restarted <- character(0)
  if(isTRUE(continue)){
    restarted <- getClonesToRestart(dir, id, all_ids,
      function(cid) clone_options(clones[match(cid, all_ids), ]), xml_md5)
    # delete the checkpoint and the BEAST files of these clones, so they start again
    cp_dir <- getCheckpointDir(dir, id)
    clone_files <- getCloneFiles(dir, id, all_ids)
    for(cid in restarted){
      unlink(file.path(cp_dir, paste0("clone_", cid, c("", ".new"))), recursive=TRUE)
      unlink(file.path(dir, clone_files[[cid]]))
    }
  }

  # In parallel runs, output of workers goes to log files, so have main
  # process show the nucleotide frequencies that the XML will use, if the
  # template uses them.
  describe <- NULL
  uses_freqs <- tryCatch(any(grepl("\\$\\{EMP_EQFREQS\\}", 
    readLines(template, warn=FALSE))), error=function(e)FALSE, warning=function(w)FALSE)
  if(uses_freqs){
    describe <- function(cs){
      freqs <- getEmpNucleotideFreqs(cs$data[[1]]@data$sequence)
      paste0("nucleotide frequencies A C G T: ", paste(freqs, collapse=" "))
    }
  }

  clones <- runTimeTreesIterate(clones, work, dir=dir, id=id, iterations=iterations, 
    ess_cutoff=ess_cutoff, ignore=ignore, continue=continue, nproc=nproc, 
    quiet=quiet, skip_errors=skip_errors, describe=describe, beast_options=clone_options)

  # add in final readBEAST with full posterior options if specified
  if(!"none" %in% posterior && nrow(clones) > 0){
    clones <- readBEAST(clones=clones, dir=dir, id=id, beast=beast, burnin=burnin, 
      trait=trait, quiet=quiet, nproc=nproc, posterior=posterior, asr=asr, 
      low_ram=low_ram, trim_ids=trim_ids)
    clones$below_ESS <- vapply(clones$parameters, countParametersBelowESS, integer(1), ignore=ignore, ess_cutoff=ess_cutoff)
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
#' @param    mcmc_length  Number of MCMC steps
#' @param    time         Name of sample time column  
#' @param    log_every    Frequency of states logged. "auto" will divide
#'                        mcmc_length by log_target         
#' @param    burnin       Burnin percent (default 10)                 
#' @param    trait        Trait column to be used. The value "?" is treated as ambiguous
#'                        (any type) and is not counted as a type.
#' @param    resume_clones  Clones to resume for \code{mcmc_length} more steps            
#' @param    include_germline Include germline sequence in analysis?     
#' @param    seq          Sequence column in data      
#' @param    germline_range   Possible date range of germline tip     
#' @param    java         Use the -java flag for BEAST run      
#' @param    seed         Use specified seeed for the -seed option for BEAST   
#' @param    log_target   Target number of samples from MCMC chain         
#' @param    tree_states  Use \code{states} vector for starting tree
#' @param    nproc      Number of cores for parallelization. At most 1 core/tree can be used.
#' @param    quiet      amount of rubbish to print to console
#' @param    rm_temp    remove temporary files (default=TRUE)
#' @param    trees      optional list of starting trees, either phylo objects or newick strings
#' @param    germline_trait_value     trait value for germline, default '?' for ambiguous
#' @param    trait_list   optional vector of all types of the trait. If NULL, it is made from the values in \code{clones}.
#' @param    ...        Additional arguments passed to tree building programs
#'
#' @return   A tibble with a column of \code{phylo} objects and \code{parameters} column
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
        tree_states=FALSE, trees=NULL, germline_trait_value='?', trait_list=NULL, ...) {

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
                            germline_trait_value=germline_trait_value,
                            trait_list=trait_list,
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
  clones$ESS_100 <- sapply(clones$parameters, function(x)sum(x$ESS < 100, na.rm=TRUE))
  clones$ESS_200 <- sapply(clones$parameters, function(x)sum(x$ESS < 200, na.rm=TRUE))

  clones
}


#' Read in a directory from a BEAST run. Runs treeannotator and loganalyser.
#' 
#' @param    data     a list of \code{airrClone} objects
#' @param    template   XML template
#' @param    beast      location of beast binary directory (beast/bin)
#' @param    dir        directory where temporary files will be placed.
#' @param    id         unique identifer for this analysis
#' @param    mcmc_length  Number of MCMC steps
#' @param    time         Name of sample time column 
#' @param    log_every    Frequency of states logged. \code{auto} will divide mcmc_length by log_target         
#' @param    burnin       Burnin percent (default 10)                 
#' @param    trait        Trait column used         
#' @param    asr          Log ancestral sequences?
#' @param    posterior    Read un full distribution of parameters and trees? Can be "none" to just have
#'  summary objects, "all" to have parameters, trees, and trees_with_traits, or a vector with the desired
#'  combination of "parameters", "trees_with_traits", and "trees".
#' @param    resume_clones  Clones to resume for \code{mcmc_length} more steps            
#' @param    include_germline Include germline in analysis?     
#' @param    start_date       Starting date of time tree if desired
#' @param    max_start_date   Maximum starting date of time tree if desired
#' @param    germline_range   Possible date range of germline tip     
#' @param    java         Use the -java flag for BEAST run      
#' @param    seed         Use specified seeed for the -seed option for BEAST   
#' @param    log_target   Target number of samples over \code{mcmc_length}         
#' @param    tree_states  Use \code{states} vector for starting tree
#' @param    nproc      Number of cores for parallelization. Uses at most 1 core per tree.
#' @param    quiet      Amount of rubbish to print to console
#' @param    low_ram    run with less memory (slightly slower)  
#' @param    trees                    optional list of starting trees, either phylo objects or newick strings
#' @param    start_edge_length        edge length to use for all branches in starting tree 
#' @param    germline_trait_value     trait value for germline, default '?' for ambiguous
#' @param    trait_list               optional vector of all types of the trait. If NULL, it is made from the values in \code{data}.
#' @param    xml_check                optional function called with the paths of the XML files after they are
#'                                    written and before BEAST starts. If it stops, BEAST is not started. 
#'                                    getTimeTreesIterate uses this to get the XML of a clone without running BEAST.
#' @param    ...      Additional arguments for XML writing functions
#'
#' @return   The input clones tibble with an additional column for the bootstrap replicate trees.
#'  
#' @seealso \link{getTimeTrees}
#' @export
buildBeast <- function(data, beast, time, template, dir, id, mcmc_length = 1000000, 
                   resume_clones=NULL, trait=NULL, asr=FALSE,
                   posterior=c("none","all","parameters","trees_with_traits","trees"),
                   log_every="auto",include_germline = TRUE, nproc = 1, quiet=0, 
                   burnin=10, low_ram=TRUE, germline_range=c(-10000,10000), java=TRUE, 
                   seed=NULL, log_target=10000, trees=NULL, tree_states=FALSE, 
                   start_edge_length=100, start_date=NULL, max_start_date=NULL,
                   germline_trait_value='?', trait_list=NULL, xml_check=NULL, ...) {

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

  # made before the clones are cut to resume_clones, so that all clones in
  # an analysis have the same types
  if(is.null(trait)){
    trait_list <- NULL
  } else if(is.null(trait_list)){
    trait_list <- makeTraitList(lapply(data, function(x)x@data[[trait]]))
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
      germline_trait_value=germline_trait_value,
      ...)

  xml_filepath <- xml_filepath[!is.na(xml_filepath)]

  # lets the caller stop before BEAST runs, for example if the XML changed
  if(!is.null(xml_check)){
    xml_check(xml_filepath)
  }

  # Run BEAST on each tree in parallel
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
  trait=trait, quiet=quiet, nproc=nproc, posterior=posterior, asr=asr, 
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
#' @param    germline_value           trait value for germline, default '?' for ambiguous
#'
#' @return   String of XML of the trait or traitSet
#'  
create_traitset <- function(clone, trait_name, column, id, trait_data_type=NULL, 
  isSet=FALSE, include_germline_as_tip=FALSE, germline_value='?') {

  all_traits <- paste(clone@data$sequence_id, clone@data[[column]], 
    collapse=",\n", sep="=")
  if (include_germline_as_tip) {
    all_traits <- paste(all_traits, paste0('Germline','=', germline_value), sep=",\n")
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
    ntips <- length(tree$tip.label)
    if(tree_states){
      tree$node.label <- (ntips + 1): (ntips + 1 + tree$Nnode - 1)
      tree$edge.length <- rep(start_edge_length, length(tree$edge.length)) # set all edge lengths to 10
    }else{
      tree$node.label <- NULL
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
#' @param    germline_trait_value     trait value for germline, default '?' for ambiguous
#' @param    root_trait               trait value of the root, if a fixed root state desired
#' @param    ...                      additional arguments for XML writing functions
#'
#' @return   File path of the written XML file
#'  
write_clone_to_xml <- function(clone, file, id, time=NULL, trait=NULL, 
  trait_data_type=NULL, template=NULL, mcmc_length=1000000, log_every=1000, replacements=NULL, 
  include_germline_as_root=FALSE, include_germline_as_tip=FALSE, 
  germline_range=c(-10000,10000), tree=NULL, trait_list=NULL, log_every_trait=10, tree_states=FALSE,
  start_edge_length=100, start_date=NULL, max_start_date=NULL, germline_trait_value='?', root_trait=NULL, ...) {
  
  kwargs <- list(...)

  # read in a template file
  if (is.null(template)) {
    # template <- system.file("extdata", "template.xml", package = "scoper")
    template <- "template.xml"
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
      trait_data_type, isSet=TRUE, include_germline_as_tip=include_germline_as_tip, germline_value=germline_trait_value)
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
      mrca_priors <- paste0(create_MRCA_prior_observed(clone, id), "\n",
        create_MRCA_prior_germline(clone, id, germline_range), "\n")
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
      if(any(grepl("GermlineRootTree", xml))) {
        warning("GermlineRootTree operator found in template, skipping TipDatesRandomWalker operator addition")
        operators <- ""
      }
      else {
      operators <-  
      paste0('<operator id="TipDatesRandomWalker.01" windowSize="1" spec="beast.base.evolution.operator.TipDatesRandomWalker" taxonset="@germSet" tree="@Tree.t:',
        id, '_', clone@clone,'" weight="1.0"/>\n')
      }
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
    freqs <- getEmpNucleotideFreqs(clone@data$sequence)
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
        # replace the value of the nodeTypes parameter with the starting traits
        node_types_index <- grep("id='nodeTypes'", xml)
        if (length(node_types_index) > 0) {
          xml[node_types_index] <- 
            gsub('value="[^"]*"', 
                 paste0('value="', paste(starting_traits, collapse=" "), '"'), 
                 xml[node_types_index])
        } else {
          stop("Could not find nodeTypes parameter in the template file")
        }
      }
    } else {
      stop("Could not find <init> tag in the template file")
    }
  }

  if (!is.null(root_trait)) {
    if (any(grepl("\\$\\{ROOT_TYPE_PROBABILITIES\\}", xml))) {
      if (is.numeric(root_trait)) {
        root_trait_index <- root_trait + 1 # add 1 because BEAST uses 0-based indexing for traits
      } else {
        root_trait_index <- match(root_trait, trait_list)
        if (is.na(root_trait_index)) {
          stop(paste0("root_trait ", root_trait, " not found in trait_list"))
        }
      }
      root_trait_probabilities <- rep(0, length(trait_list))
      root_trait_probabilities[root_trait_index] <- 1
      root_trait_probabilities_string <- paste0(root_trait_probabilities, collapse=" ")
      root_trait_probabilities_string <- paste0('<rootFrequencies id="rootfreqs.s:newTrait" spec="Frequencies">
                        <parameter id="rootFrequencies.s:newTrait" spec="parameter.RealParameter" dimension="', length(trait_list), '" name="frequencies">', root_trait_probabilities_string, '</parameter>
                    </rootFrequencies>')
      xml <- gsub("\\$\\{ROOT_TYPE_PROBABILITIES\\}", root_trait_probabilities_string, xml)
    } else {
      warning("root_trait argument provided but ${ROOT_TYPE_PROBABILITIES} placeholder not found in template file")
    }
  } else {
    if (any(grepl("\\$\\{ROOT_TYPE_PROBABILITIES\\}", xml))) {
      xml <- gsub("\\$\\{ROOT_TYPE_PROBABILITIES\\}", "", xml)
    }
  }
  
  curr_replace_name = "NODES_TYPE_INIT"
  if (any(grepl("\\$\\{NODES_TYPE_INIT\\}", xml))) {
    if (!curr_replace_name %in% names(kwargs)) {
      nodes_type_init = 0
    } else {
      nodes_type_init = kwargs[[curr_replace_name]]
    }
    xml <- gsub(paste0("\\$\\{", curr_replace_name, "\\}"), nodes_type_init, xml)
  }

  matches <- unlist(regmatches(xml, gregexpr("\\$\\{([^}]+)\\}", xml)))
  template_variables <- unique(sub("\\$\\{([^}]+)\\}", "\\1", matches))
  # print(paste("Template variables found:", paste(template_variables, collapse=", ")))

  for (var in template_variables) {
    if (!var %in% names(kwargs)) {
      msg <- paste("Variable", var, "not provided but is required by the template.")
      hint <- .templateVarLookup(xml, var)
      if (!is.null(hint) && nchar(hint) > 0) {
        msg <- paste0(msg, " (", hint, ")")
      } else {
        msg <- paste0(msg, " Run describeTemplate() on this template file for a full list of what it needs.")
      }
      stop(msg)
    }
    if (var %in% replacements) {
      next # skip replacements, they are handled later
    }

    value <- kwargs[[var]]

    # replace the ${VAR} placeholder with the value from kwargs
    xml <- gsub(paste0("\\$\\{", var, "\\}"), value, xml)
  }

  # strip out any template info comments from the xml as they are not useful in the final xml file, but do it before replacements since those are often large and we want to keep string operations involving them to a minimum
  xml <- .stripTemplateInfo(xml)

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

#' Empirical frequencies of the nucleotides A, C, G and T in sequences
#'
#' @param    seqs   vector of sequences
#'
#' @return   Table of the frequencies of A, C, G and T, in this order
#' @noRd
getEmpNucleotideFreqs <- function(seqs){
  nts <- c("A", "C", "G", "T")
  splits <- strsplit(seqs, split="")
  counts <- table(unlist(splits))[nts]
  counts/sum(counts)
}

#' Make the sorted list of types for a trait from the values in all clones
#' 
#' In BEAST XMLs, "?" is ambiguous, so should not be counted as a type itself.
#'
#' @param    values   list of vectors of trait values, one vector per clone
#'
#' @return   Sorted vector of the unique trait values, without "?"
#' @noRd
makeTraitList <- function(values){
  values <- setdiff(sort(unique(unlist(values))), "?")
  if(length(values) == 0){
    stop("No trait values found other than '?'. At least one type is required.")
  }
  values
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
#' @param    germline_trait_value     trait value for germline, default '?' for ambiguous
#' @param    ...                      additional arguments for XML writing functions
#'
#' @return   File paths of the written XML files
#'  
write_clones_to_xmls <- function(data, id, trees=NULL, time=NULL, trait=NULL, template=NULL, 
  outfile=NULL, replacements=NULL, trait_list=NULL, 
  mcmc_length=1000000, log_every=1000, include_germline_as_root=FALSE, 
  include_germline_as_tip=FALSE, germline_range=c(-10000,10000), 
  tree_states=FALSE, start_edge_length=100, start_date=NULL, max_start_date=NULL, germline_trait_value='?',
  ...) {

  kwargs <- list(...)

  # iterate over the clones to first create trait data type if trait exists
  if (!is.null(trait)) {
    if (is.null(trait_list)) {
      # get all the possible values of the trait
      trait_list <- makeTraitList(lapply(data, function(x)x@data[[trait]]))
    } else {
      # "?" is never a type, even if the caller put it in the list
      if ("?" %in% trait_list) {
        warning("Trait list contains '?', which is not a valid type. Removing it from the list.")
        trait_list <- setdiff(trait_list, "?")
      }
    }
    codeMap <- paste(trait_list, 0:(length(trait_list)-1), collapse=",\n", sep="=")
    codeMap <- paste0(codeMap, ",\n? = ", paste0(0:(length(trait_list)-1), collapse=" "))
    trait_data_type <- paste0('<userDataType id="traitDataType.newTrait" spec="beast.base.evolution.datatype.UserDataType" codeMap="',
     codeMap, '" codelength="-1" states="', length(trait_list), '"/>')
  }
  xmls <- c()
  for (i in 1:length(data)){
    #if ("trees" %in% names(kwargs)) {
    if(!is.null(trees)){
      # if trees are provided, use them
      tree <- trees[[i]]
    } else {
      tree <- NULL
    }
    xmls <- c(xmls, write_clone_to_xml(data[[i]], 
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
                     germline_trait_value=germline_trait_value,
                     ...))
  }
  return(xmls)
}


#' Reads in a BEAST output directory
#' 
#' \code{readBEAST} Reads in data from BEAST output directory
#' @param clones     either a tibble (getTrees) or list of \code{airrClone} objects
#' @param beast      location of beast binary directory (beast/bin)
#' @param dir        directory where BEAST output files have been placed.
#' @param id         unique identifer for this analysis
#' @param trait      Trait column used         
#' @param asr        Log ancestral sequences?
#' @param posterior    Read un full distribution of parameters and trees? Can be "none" to just have
#'  summary objects, "all" to have parameters, trees, and trees_with_traits, or a vector with the desired
#'  combination of "parameters", "trees_with_traits", and "trees".
#' @param nproc      Number of cores for parallelization. Uses at most 1 core per tree.
#' @param quiet      amount of rubbish to print to console
#' @param burnin         percent of initial tree samples to discard (default 10)
#' @param low_ram        run with less memory (slightly slower)
#' @param trim_ids    remove last _ group from tips?       
#'
#' @return   
#' If data is a tibble, then the input clones tibble with additional columns for 
#' trees and parameter estimates given the specified burnin. If input is just a 
#' list of airrClone objects, it will return the corresponding list of trees
#' given the burnin
#'  
#' @export
readBEAST <- function(clones, dir, id, beast, burnin=10, trait=NULL, nproc = 1, 
  quiet=0, posterior=c("none","all","parameters","trees_with_traits","trees"), 
  asr=FALSE, low_ram=TRUE, trim_ids=FALSE) {

  if(!"list" %in% class(clones) && "data" %in% names(clones)){
    data <- clones$data
  }else if("list" %in% class(clones)){
    data <- clones
  }else{
    stop("Input data type not supported")
  }

  if(setequal(posterior, c("all","none","parameters","trees_with_traits","trees"))){
    # default
    posterior = c("parameters","trees_with_traits","trees")
  }else if("all" %in% posterior){
    posterior = c("parameters","trees_with_traits","trees")
  }else if("none" %in% posterior){
    posterior = "none"
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
  trees <- list()
  for(i in 1:length(data)){
    if(is.null(trait)){
      treefile <- file.path(dir, paste0(id,"_",data[[i]]@clone, ".tree"))
    }else{
      treefile <- file.path(dir, paste0(id,"_",data[[i]]@clone, "_tree_with_trait.tree"))
    }
    logfile <- file.path(dir, paste0(id,"_",data[[i]]@clone, ".log"))
    logoutfile <- file.path(dir, paste0(id, "_", data[[i]]@clone,"_log.tsv"))

    beastobj <- treeio::read.beast(treefile)
    if("error" %in% class(beastobj)){
      stop(paste("Couldn't read in ",treefile))
    }

    # add the usual bells and whistles for Dowser tree objects
    beastobj@phylo$name <- data[[i]]@clone
    beastobj@phylo$tree_method <- paste("beast")
    beastobj@phylo$edge_type <- "time"
    nnodes <- length(unique(c(beastobj@phylo$edge[,1],beastobj@phylo$edge[,2])))
    beastobj@phylo$nodes <- lapply(1:nnodes, function(x)list(sequence=NA))

    # add in node sequences. Currently only adds tips.
    for(j in 1:length(beastobj@phylo$nodes)){
      if(j <= length(beastobj@phylo$tip.label)){
        label <- beastobj@phylo$tip.label[j]
        if(label != "Germline"){
          beastobj@phylo$nodes[[j]]$sequence <- dplyr::filter(data[[i]]@data, 
            !!rlang::sym("sequence_id") == label)[[data[[i]]@phylo_seq]]
        }else{
          beastobj@phylo$nodes[[j]]$sequence <- switch(data[[i]]@phylo_seq,
            "sequence" = data[[i]]@germline,
            "lsequence" = data[[i]]@lgermline,
            "hlsequence" = data[[i]]@hlgermline)
        }
      }
    }

    l <- readLines(logoutfile)
    # just in case there are extra lines at the top of the log file, find the line where the parameter table starts
    item_lines <- grep("item", l)
    if (length(item_lines) == 1) {
      log <- read.table(text=l[(item_lines):(length(l)-1)], header=TRUE)
    } else {
      log <- read.table(text=l[4:(length(l)-1)], header=TRUE)
    }

    # add parameter summary
    beastobj@info$parameters <- log
    if("trees" %in% posterior){ 
      treesfile <- file.path(dir, paste0(id,"_",data[[i]]@clone, ".trees"))      
      l <- readLines(treesfile, warn=FALSE)
      if(!grepl("End;",l[length(l)])){
        l[length(l) + 1] <- "End;"
        warning("Adding End; to ",treesfile)
        #make new file to avoid overwriting
        treesfile <- file.path(dir, paste0(data[[i]]@clone, "_end.trees"))
        writeLines(l, con=treesfile)
      }
      # always gives an incomplete final line warning, but not 
      # an issue if we end on End;
      phylos <- suppressWarnings(treeio::read.beast(treesfile))
      burn <- floor(length(phylos)*burnin/100)
      phylos <- phylos[(burn+1):length(phylos)]
      if(trim_ids){
        phylos <- lapply(phylos, function(y){
        y@phylo$tip.label <- sapply(strsplit(y@phylo$tip.label,"_"), function(x)
          paste0(x[1:(length(x)-1)], collapse="_"))
        y
      })}
      
      beastobj@info$trees_posterior <- phylos
    }
    if("trees_with_traits" %in% posterior){ 
      if(is.null(trait)){
        stop("trait column must be specified when trees_with_traits %in% posterior")
      }else{
        treesfile <- file.path(dir, paste0(id,"_",data[[i]]@clone, "_tree_with_trait.trees"))
      }
      l <- readLines(treesfile, warn=FALSE)
      if(!grepl("End;",l[length(l)])){
        l[length(l) + 1] <- "End;"
        warning("Adding End; to ",treesfile)
        #make new file to avoid overwriting
        treesfile <- file.path(dir, paste0(data[[i]]@clone, "_end.trees"))
        writeLines(l, con=treesfile)
      }
      phylos <- suppressWarnings(treeio::read.beast(treesfile))
      burn <- floor(length(phylos)*burnin/100)
      phylos <- phylos[(burn+1):length(phylos)]
      if(trim_ids){
        phylos <- lapply(phylos, function(y){
        y@phylo$tip.label <- sapply(strsplit(y@phylo$tip.label,"_"), function(x)
          paste0(x[1:(length(x)-1)], collapse="_"))
        y
      })}
      beastobj@info$trees_with_traits_posterior <- phylos
    }
    if("parameters" %in% posterior){
      l <- read.table(logfile, header=TRUE)
      burn <- floor(nrow(l)*burnin/100)
      l <- l[(burn+1):nrow(l),]
      beastobj@info$parameters_posterior <- l
    #  beastobj@info$parameters_posterior <- tidyr::gather(l,
    #   "parameter", "value", -(!!rlang::sym("Sample")))
    }
    beastobj@info$name <- data[[i]]@clone
    trees[[i]] <- beastobj
  }

  # set divergence of each node
  for(i in 1:length(trees)){
    trees[[i]]@phylo <- setNodeDivergences(trees[[i]]@phylo)
  }

  if(quiet < 1)print("Ran readBEAST")

  if(!"list" %in% class(clones) && "data" %in% names(clones)){
    clones$trees <- trees
    clones$parameters <- lapply(trees, function(x)x@info$parameters)
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
#' @param  object    treedata object with parameters_posterior and trees_posterior
#' @param  bins      number of bins for plotting
#' @param  youngest  timepoint of the most recently tip sampled (if 0, backward time used)
#' @param  clone_id  name of the clone being analyzed (if desired)
#' @param  max_height max height to use (min, median, mean, max)
#' @param  exclude_germline exclude germline from skyline plot? (For TyCHE germline-rooted trees)
#' @return   Bayesian Skyline values for given clone
#'
#' @export
makeSkyline <- function(object, bins=100, youngest=0, 
    clone_id=NULL, max_height=c("min","median","mean","max"), exclude_germline=TRUE){
    
    if(length(max_height) > 1){
      max_height <- max_height[1]
    }
    params <- object@info$parameters_posterior
    phylos <- object@info$trees_posterior
    if(is.null(params) || is.null(phylos)){
      stop(paste("Parameters and trees posterior not found.\n",
        "Try reading in object with readBEAST(posterior='all')"))
    }
    phylos <- lapply(phylos, function(x)x@phylo)
    params <- tidyr::gather(params, "parameter", "value", -(!!rlang::sym("Sample")))

    if(!"bPopSizes.1" %in% unique(params$parameter)){
        stop(paste("log file doesn't have pop sizes.",
            "Was it run with skyline tree_prior='coalescent_skyline'?"))
    }

    samples <- sort(unique(params$Sample))
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
        stop(paste(object@phylo$name, "found popsizes < 0, can't continue"))
    }
    if(sum(groups$value < 0) > 0){
        stop(paste(object@phylo$name, "found groupsizes < 0, can't continue"))
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
      if (exclude_germline){
        germline_name <- ifelse("Germline" %in% tr$tip.label, "Germline", 
          ifelse("germline" %in% tr$tip.label, "germline", 
            ifelse("GL" %in% tr$tip.label, "GL", NA)))
        if(is.na(germline_name)) {
          warning(paste( 
            "Couldn't find germline tip in tree, proceeding without dropping germline"))
        } else {
          tr <- ape::drop.tip(tr, germline_name)
        }
      }
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
        warning(paste(object@phylo$name, "Removing",length(indistinct),
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

    return(as.data.frame(skyplot))
}

#' Make data frames for Bayesian skyline plots
#' 
#' \code{makeSkylines} 
#' @param  clones    clone tibble
#' @param  time      name of time column
#' @param  bins      number of bins for plotting
#' @param  verbose   if 1, print name of clones
#' @param  forward   plot in forward or (FALSE) backward time?
#' @param  nproc     processors for parallelization (by clone)
#' @param  max_height max height to use (min, median, mean, max)
#' @param  exclude_germline exclude germline from skyline plot? (For TyCHE GRTs)
#' @return   Bayesian Skyline values for given clone
#' @details Clones must contain treedata objects with parameters_posterior and
#' trees_posterior. See \code{readBEAST} or \code{getTimeTrees} with posterior="all"
#' @export
getSkylines <- function(clones, time, bins=100, verbose=0, forward=TRUE,
    nproc=1, max_height=c("min","median","mean","max"), exclude_germline=TRUE){

    if(is.null(clones$trees[[1]]@info$parameters_posterior) || 
      is.null(clones$trees[[1]]@info$trees_posterior)){
      stop(paste("Parameters and trees posterior not found.\n",
        "Try reading in object with readBEAST(posterior='all')"))
    }
    if(forward){
        youngest <- sapply(clones$data, function(x)
            max(as.numeric(x@data[[time]])))
    }else{
        youngest <- rep(0, length=nrow(clones))
    }

    skylines <- parallel::mclapply(1:nrow(clones), function(x){
        if(verbose != 0){
            print(paste(clones$clone_id[x], youngest[x]))
        }
        tryCatch(makeSkyline(clones$trees[[x]],
            youngest=youngest[x], bins=bins, 
            clone_id=clones$clone_id[x], max_height=max_height, 
            exclude_germline=exclude_germline),
            error=function(e)e)
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

#' Add height and length columns to a tree@data data frame
#' @param  tree    a treedata objects from read.beast
#' @export
getHeightsAndLengths = function(tree){
  phylo <- tree@phylo
  uca <- ape::getMRCA(phylo, tip=phylo$tip.label)
  dists <- ape::dist.nodes(phylo)
  max_height <- max(dists[uca,])
  tree@data$height <- NA
  tree@data$length <- NA
  for(i in 1:nrow(tree@data)){
    node <- tree@data$node[i]
    tree@data$height[i] <- max_height - dists[node,uca]
    if(node == uca){
      tree@data$length[i] <- 0
    }else{
      tree@data$length[i] <- phylo$edge.length[which(phylo$edge[,2] == node)]
    }
  }
  return(tree)
}


#' Recurse up to tree to find most recent node with different state, or the root
#' 
#' \code{getDiffPoint} 
#' @param  tree    a treedata object, from getTimeTrees
#' @param  targetnode    current node
#' @param  trait   column name of the trait of interest
#' @param  height  which height value to return
#' @param verbose  print out run info
#' @param eo_adjust adjust heights using expectOccupancies (recommended, requires eo_type)
#' @param eo_type  if eo_adjust, trait value described by expectedOccupancies (typically state 1 of 2)
#' 
#' @export
getDiffPoint = function(tree, targetnode, trait, height="height", verbose=FALSE,
  eo_adjust=FALSE, eo_type=NULL){
  type <- dplyr::filter(tree@data, !!rlang::sym("node")==targetnode)[[trait]]
  edge <- tree@phylo$edge[tree@phylo$edge[,2] == targetnode,]
  if(verbose){
    print(paste(targetnode,type))
    print(edge)
  }
  if(length(edge) == 0){
    return(dplyr::tibble(parent_node=targetnode, parent_type=type, 
      parent_height=as.numeric(filter(tree@data, !!rlang::sym("node")==targetnode)[[height]]),
      child_node=targetnode, child_type=type, 
      child_height=as.numeric(filter(tree@data, !!rlang::sym("node")==targetnode)[[height]]),
      root=1
      ))
  }
  if(!is.null(nrow(edge))){
    stop("weird")
  }
  parent <- as.character(edge[1])
  parent_type <- dplyr::filter(tree@data, !!rlang::sym("node")==parent)[[trait]]
  parent_height <- as.numeric(dplyr::filter(tree@data, !!rlang::sym("node")==parent)[[height]])
  child_height <- as.numeric(dplyr::filter(tree@data, !!rlang::sym("node")==targetnode)[[height]])
  if(parent_type == type){
    return(getDiffPoint(tree=tree, targetnode=parent, trait=trait, 
      height=height, verbose=verbose, eo_adjust=eo_adjust, eo_type=eo_type))
  }else{
    if(eo_adjust){
      child <- dplyr::filter(tree@data, !!rlang::sym("node")==targetnode)
      if(type == eo_type){
        parent_height <- parent_height - 
          (1-as.numeric(child$expectedOccupancies)) * 
          as.numeric(child$length)
        child_height <- child_height + 
          (as.numeric(child$expectedOccupancies)) * 
          as.numeric(child$length)
      }else{ # if parent is EO type (i.e. tip is not, assumes only 2 types)
        parent_height <- parent_height - 
          as.numeric(child$expectedOccupancies) * 
          as.numeric(child$length)
        child_height <- child_height + 
          (1-as.numeric(child$expectedOccupancies)) * 
          as.numeric(child$length)
      }
    }
    return(dplyr::tibble(parent_node=parent, 
      parent_type=parent_type, parent_height=parent_height,
      child_type=type, child_height=child_height, child_node=targetnode,
      root=0))
  }
}

#' For each tree, recurse up to tree to find most recent 
#' node with a different state, or the root
#' 
#' \code{getDiffPoints} 
#' @param  data    a tibble containing trees column from getTimeTrees
#' @param  trait   column name of the trait of interest
#' @param  height  which height value to return
#' @param  verbose print out info during run
#' @param  tip_traits vector of other traits to include for each tip.
#'           Must have been also specified as traits in formatClones.
#' @param eo_adjust adjust heights using expectOccupancies. Recommended if EO model used,
#'         requires eo_type to specify the type whose occupancy is tracked)
#' @param eo_type  if eo_adjust, trait value described by expectedOccupancies (typically state 1 of 2)
#' @param full_posterior Computer statistics using the full posterior distribution of trees?
#' @param summarize if full_posterior=TRUE, summarize results or return the full table?
#' @param nproc Number of cores to use (parallelizes by row of input data object)
#' @details 
#' Returns a data frame where each row is a tip in each tree
#' clone_id = clone id for that tree
#' tip = tip name
#' tip_tip = trait value for that tip
#' tip_height = height value for that tip
#' parent_node = most recent ancestor node with different trait value, or root
#' parent_type = type of parent_node, will be "root" if at root node
#' parent_height = height of parent_node 
#' child_node = most recent ancestor node before parent_node, or root
#' child_type = type of child_node, will be "root" if at root node
#' child_height = height of child_node 
#' root = 1 if at root node, 0 otherwise
#' #' <other columns> copied over from airrClone object for each tip
#' 
#' if summarize=TRUE, parent_node and child_node will not be available
#' instead:
#' tip_type = type of tip
#' parent_type = type of parent at different node
#' samples = number of samples using this tip, tip_type, parent_type
#' tip_height = height of tip           
#' parent_height_mean = mean of parent_height
#' parent_height_95HPDlo = lower 95% HPD interval of parent_height
#' parent_height_95HPDup = upper 95% HPD interval of parent_height
#' child_height_mean = mean of parent_height     
#' child_height_95HPDlo = lower 95% HPD interval of parent_height  
#' child_height_95HPDup = upper 95% HPD interval of parent_height  
#' root_freq = frequency at which root == 1            
#' 
#' @examples
#' 
#' \dontrun{dp = getDiffPoints(data, "location", verbose=TRUE, 
#'  eo_adjust=TRUE, eo_type="germinal_center")}
#' 
#' @export
getDiffPoints = function(data, trait, height="height", verbose=FALSE,
  tip_traits=NULL, eo_adjust=FALSE, eo_type=NULL, full_posterior=FALSE,
  summarize=TRUE, nproc=1){
  
  results_list <- parallel::mclapply(1:nrow(data),function(x){
    if(verbose)print(data$clone_id[x])
    trees <- data$trees[[x]]
    if(full_posterior){
      if(!"trees_with_traits_posterior" %in% names(trees@info)){
        stop(paste("Tree posterior distribution not found.\nEither first run readBEAST using",
          "posterior='trees_with_traits' (or posterior='all') or set full_posterior=FALSE to use a single tree"))
      }
      trees <- trees@info$trees_with_traits_posterior
    }else{
      trees <- list(trees)
    }
    diffpoints <- dplyr::tibble()
    treecounter <- 1
    for(tree in trees){
      if(!"height" %in% names(tree@data) || !"length" %in% names(tree@data)){
        tree <- getHeightsAndLengths(tree)
      }
      for(l in tree@phylo$tip.label){
        #if(verbose){
        #  print(paste(treecounter,l))
        #}
        d <- dplyr::filter(tree@data, !!rlang::sym("node") == 
          which(tree@phylo$tip.label == l))
        df <- getDiffPoint(tree, targetnode=as.character(which(tree@phylo$tip.label == l)), 
          trait=trait, height=height, verbose=FALSE, eo_adjust=eo_adjust,
          eo_type=eo_type)
        temp <- dplyr::tibble(clone_id=data$clone_id[x], tip=l, 
          tip_type=d[[trait]], tip_height=d[[height]], treecounter=treecounter)
        diffpoints <- dplyr::bind_rows(diffpoints, dplyr::bind_cols(temp, df))
        treecounter <- treecounter + 1
      }
    }
    diffpoints$parent_height <- as.numeric(diffpoints$parent_height)
    diffpoints$child_height <- as.numeric(diffpoints$child_height)
    diffpoints$tip_height <- as.numeric(diffpoints$tip_height)

    getHPD = function(x){
      if(length(x) == 1){
        warning("Singular value in HPD calculation")
        return(c(x,x))
      }else{
        return(coda::HPDinterval(coda::as.mcmc(x), prob = 0.95))
      }
    }
    
    if(summarize & length(trees) > 1){
      diffpoints <- diffpoints %>%
        group_by(!!rlang::sym("tip"), !!rlang::sym("tip_type"), !!rlang::sym("parent_type")) %>%
        summarize(
          samples = n(),
          tip_type = unique(!!rlang::sym("tip_type")),
          tip_height = mean(!!rlang::sym("tip_height")),
          parent_height_mean = mean(!!rlang::sym("parent_height")),
          parent_height_95HPDlo = getHPD(!!rlang::sym("parent_height"))[1],
          parent_height_95HPDup = getHPD(!!rlang::sym("parent_height"))[2],
          child_height_mean = mean(!!rlang::sym("child_height")),
          child_height_95HPDlo = getHPD(!!rlang::sym("child_height"))[1],
          child_height_95HPDup = getHPD(!!rlang::sym("child_height"))[2],
          root_freq = mean(!!rlang::sym("root")),
          .groups = "drop_last"
          )
    }

    # copy over trait info for each tip
    if(!is.null(tip_traits)){
      for(tr in tip_traits){
        if(!tr %in% names(data$data[[x]]@data)){
          stop(paste(tr, "not found in airrClone object"))
        }
        diffpoints[[tr]] = NA
        for(tindex in 1:nrow(diffpoints)){
          m <- match(diffpoints$tip[tindex], data$data[[x]]@data$sequence_id)
          if(is.na(m)){
            if(diffpoints$tip[tindex] == "Germline"){
              diffpoints[[tr]][tindex] <- NA
              next
            }
            stop(paste(diffpoints$tip[tindex],"not found in airrClone object"))
          }
          diffpoints[[tr]][tindex] <- data$data[[x]]@data[[tr]][m]
        }
      }
    }
    diffpoints$clone_id <- data$clone_id[x]
    diffpoints
  }, mc.cores=nproc)
  results <- tryCatch(dplyr::bind_rows(results_list), 
    error=function(e){
      saveRDS(results_list, "getDiffPoints_error.rds")
      stop("Error caught in getDiffPoints, see getDiffPoints_error.rds")
    })
  return(results)
}


#' Parse a dowser-template-info comment block out of a template's raw lines
#'
#' @param    xml   character vector, one element per line (as from
#'                  \code{readLines})
#' @return   a list with \code{$vars} (named list, each element a list with
#'           \code{tier} and \code{description}), \code{$structs} (named
#'           list of free-text notes), \code{$description} (character
#'           scalar, or NULL if no \code{@description} entry was found), and
#'           \code{$found} (logical)
#' @keywords internal
.parseTemplateInfo <- function(xml) {
 
  empty <- list(vars = list(), structs = list(), description = NULL, found = FALSE)
 
  start_idx <- grep("^\\s*<!--\\s*dowser-template-info\\s*$", xml)
  if (length(start_idx) == 0) {
    return(empty)
  }
  start_idx <- start_idx[1]
 
  end_idx <- grep("^\\s*-->\\s*$", xml)
  end_idx <- end_idx[end_idx > start_idx]
  if (length(end_idx) == 0) {
    warning("Found a dowser-template-info block with no closing '-->'; ignoring it.")
    return(empty)
  }
  end_idx <- end_idx[1]
 
  if (end_idx == start_idx + 1) {
    return(list(vars = list(), structs = list(), description = NULL, found = TRUE)) # empty block
  }
  block <- xml[(start_idx + 1):(end_idx - 1)]
 
  vars <- list()
  structs <- list()
  description <- NULL
  last_entry <- NULL # list(type="var"|"struct"|"description", name=...)
 
  var_pattern         <- "^\\s*@var\\s+(\\S+)\\s*\\|\\s*([^|]+?)\\s*\\|\\s*(.*)$"
  struct_pattern      <- "^\\s*@struct\\s+(\\S+)\\s*\\|\\s*(.*)$"
  description_pattern <- "^\\s*@description\\s*\\|\\s*(.*)$"
 
  for (line in block) {
    if (grepl("^\\s*$", line)) {
      next # blank lines don't break a continuation, but add nothing either
    }
    if (grepl(var_pattern, line)) {
      m <- regmatches(line, regexec(var_pattern, line))[[1]]
      name <- m[2]
      vars[[name]] <- list(tier = trimws(m[3]), description = trimws(m[4]))
      last_entry <- list(type = "var", name = name)
    } else if (grepl(struct_pattern, line)) {
      m <- regmatches(line, regexec(struct_pattern, line))[[1]]
      name <- m[2]
      structs[[name]] <- trimws(m[3])
      last_entry <- list(type = "struct", name = name)
    } else if (grepl(description_pattern, line)) {
      m <- regmatches(line, regexec(description_pattern, line))[[1]]
      description <- trimws(m[2])
      last_entry <- list(type = "description")
    } else if (!is.null(last_entry)) {
      extra <- trimws(line)
      if (last_entry$type == "var") {
        vars[[last_entry$name]]$description <-
          paste(vars[[last_entry$name]]$description, extra)
      } else if (last_entry$type == "struct") {
        structs[[last_entry$name]] <- paste(structs[[last_entry$name]], extra)
      } else {
        description <- paste(description, extra)
      }
    }
    # a malformed line with no prior entry to attach to is silently skipped;
    # describeTemplate()'s "undocumented" reporting will surface the gap
    # for any real ${VAR} that ends up with no parsed entry.
  }
 
  list(vars = vars, structs = structs, description = description, found = TRUE)
}
 
#' Look up a single variable's documentation, for use in error messages
#'
#' @param    xml   character vector of template lines
#' @param    var   variable name (bare, no ${})
#' @return   the description string, or NULL if undocumented/no block found
#' @keywords internal
.templateVarLookup <- function(xml, var) {
  info <- .parseTemplateInfo(xml)
  if (!is.null(info$vars[[var]])) {
    info$vars[[var]]$description
  } else {
    NULL
  }
}
 
#' Built-in checklist of Dowser-recognized structural template elements
#'
#' Unlike free-form \code{${VAR}}s, these elements' behavior is defined by
#' Dowser's own substitution code (\code{write_clone_to_xml}), not by the
#' template author -- so this checklist works on any template, annotated or
#' not.
#'
#' @param    xml   character vector of template lines
#' @return   a list of check results, each a list with \code{name},
#'           \code{present} (logical), and \code{effect} (character)
#' @keywords internal
.structuralChecklist <- function(xml) {
 
  has <- function(pattern) any(grepl(pattern, xml))
 
  checks <- list(
    list(
      name = "<init> tag",
      present = has("<init") && has("</init>"),
      if_present = "Custom starting trees (the `tree` argument) are supported.",
      if_absent  = "Passing `tree=` will error (\"Could not find <init> tag in the template file\"). Fine if you never supply a starting tree -- every shipped template has this tag."
    ),
    list(
      name = "nodeTypes parameter (id='nodeTypes')",
      present = has("id='nodeTypes'"),
      if_present = "Starting trees with pre-specified type states (`tree_states=TRUE`) are supported.",
      if_absent  = "`tree_states=TRUE` with a starting tree will error (\"Could not find nodeTypes parameter in the template file\"). None of the shipped xml-templates base files declare this -- add a `<parameter id='nodeTypes' .../>` node if you need it."
    ),
    list(
      name = "GermlineRootTree",
      present = has("GermlineRootTree"),
      if_present = "Dowser will NOT add a germline TipDatesRandomWalker operator when `include_germline_as_tip=TRUE` -- this template's own GRT operators handle tip/root dating jointly (Dowser prints a warning noting the skip).",
      if_absent  = "If `include_germline_as_tip=TRUE`, Dowser will add a TipDatesRandomWalker operator via the OPERATORS placeholder, if present (see below)."
    ),
    list(
      name = "MRCA placeholder",
      present = has("\\$\\{MRCA\\}"),
      if_present = "`start_date`, `max_start_date`, and germline-tip MRCA priors (`include_germline_as_tip=TRUE`) will take effect.",
      if_absent  = "`start_date`, `max_start_date`, and germline-tip MRCA priors (`include_germline_as_tip=TRUE`) will be silently ignored if supplied -- no error is raised, no MRCA prior is added."
    ),
    list(
      name = "ROOTFREQS placeholder",
      present = has("\\$\\{ROOTFREQS\\}"),
      if_present = "`include_germline_as_root=TRUE` will take effect (root nucleotide frequencies enforced via the external rootfreqs package).",
      if_absent  = "`include_germline_as_root=TRUE` will be silently ignored -- the germline will NOT be enforced as the root state."
    ),
    list(
      name = "OPERATORS placeholder",
      present = has("\\$\\{OPERATORS\\}"),
      if_present = "A germline tip-dating operator can be added when `include_germline_as_tip=TRUE` (unless GermlineRootTree is present, see above).",
      if_absent  = "No germline tip-dating operator will ever be added, even if `include_germline_as_tip=TRUE`."
    ),
    list(
      name = "ROOT_TYPE_PROBABILITIES placeholder",
      present = has("\\$\\{ROOT_TYPE_PROBABILITIES\\}"),
      if_present = "`root_trait=` will take effect (root type pinned via a one-hot prior).",
      if_absent  = "`root_trait=` will raise a warning and be ignored -- this is the one case in this checklist that warns rather than silently doing nothing."
    ),
    list(
      name = "EMP_EQFREQS placeholder",
      present = has("\\$\\{EMP_EQFREQS\\}"),
      if_present = "Empirical per-clone nucleotide frequencies are computed from your sequence data at runtime.",
      if_absent  = "Nucleotide frequencies come from whatever is already fixed or estimated in the template; nothing is computed from your data."
    ),
    list(
      name = "DATE placeholder",
      present = has("\\$\\{DATE\\}"),
      if_present = "`time=` sample dates will be attached to tips. (If `time=` is NULL, this template will error.)",
      if_absent  = "`time=` will be silently ignored if supplied -- no dates attached, no error."
    ),
    list(
      name = "TRAIT / TRAIT_NAME placeholder",
      present = has("\\$\\{TRAIT\\}"),
      if_present = "`trait=` (with `trait_data_type=`) will be attached and reconstructed. If TRAIT_NAME is also present the trait column name is substituted directly; if only TRAIT is present, Dowser falls back to replacing tag=\"location\" (older-style templates). (If `trait=` is NULL, this template will error.)",
      if_absent  = "`trait=` will be silently ignored if supplied -- no discrete trait reconstructed, no error."
    ),
    list(
      name = "NODES placeholder",
      present = has("\\$\\{NODES\\}"),
      if_present = "Per-node parameter dimensioning for type-linked models (e.g. nodeTypes, expectedOccupancies) is filled in automatically.",
      if_absent  = "Not used by this template -- expected for non-type-linked (StrictClock/UCLD) templates."
    ),
    list(
      name = "NODES_TYPE_INIT placeholder",
      present = has("\\$\\{NODES_TYPE_INIT\\}"),
      if_present = "Initial internal-node type index at MCMC start can be set via the `NODES_TYPE_INIT` argument; defaults to 0 if omitted even when this placeholder is present.",
      if_absent  = "Not used by this template."
    )
  )
 
  checks
}
 
#' Describe an xml-templates file: what it needs and what it supports
#'
#' Reports every free-form \code{${VAR}} placeholder in a template, with
#' descriptions pulled from an embedded \code{dowser-template-info} comment
#' block if the template author provided one, plus a built-in checklist of
#' structural elements Dowser recognizes (MRCA, ROOTFREQS, OPERATORS,
#' ROOT_TYPE_PROBABILITIES, a \code{nodeTypes} parameter, and whether
#' the template uses \code{GermlineRootTree}). The structural checklist
#' works on any template, documented or not, since that behavior is defined
#' by Dowser's own code rather than by the template author -- several of
#' these are silent no-ops rather than errors when an argument is supplied
#' but the template doesn't support it, which is exactly what this function
#' is meant to catch before a long MCMC run rather than after. If the
#' template's \code{dowser-template-info} block has an overall
#' \code{@description} entry, that's printed first, before the structural
#' checklist and variable list.
#'
#' @param    template   path to an xml-templates file
#' @param    print      print a human-readable summary? Default TRUE. Set
#'                       FALSE if you only want the returned object.
#' @param    var        optional: a single variable name (bare, no
#'                       \code{${}}), free-form or Dowser-special (e.g.
#'                       \code{"EMP_EQFREQS"}). If supplied, instead of the
#'                       full dump, this looks up just that one variable and
#'                       prints/returns a focused answer: whether it's in
#'                       the template at all, and either its
#'                       \code{dowser-template-info} description (free-form)
#'                       or its built-in structural effect (Dowser-special).
#'
#' @return   invisibly: if \code{var} is NULL, a list with
#'           \code{$description} (character scalar, or NULL if the template
#'           has no \code{@description}), \code{$vars} (data.frame: name,
#'           tier, documented, description), and \code{$structural}
#'           (data.frame: element, present, effect). If \code{var} is
#'           supplied, a list with \code{$name}, \code{$kind} ("free",
#'           "structural", or "not_found"), \code{$present}, and whichever
#'           of \code{$documented}/\code{$tier}/\code{$description} (free)
#'           or \code{$effect} (structural) applies.
#'
#' @details
#' See the xml-templates variant reference at https://tyche.readthedocs.io
#' for the \code{dowser-template-info} format if you're documenting your own
#' template.
#'
#' @examples
#' \dontrun{
#' describeTemplate("TypeLinked/TypeLinkedExpectedOccupancy_EstTraitClockRates_EmpFreq.xml")
#' describeTemplate("TypeLinked/TypeLinkedExpectedOccupancy_EstTraitClockRates_EmpFreq.xml",
#'                   var = "TRAIT_RATE_MEAN_1")
#' }
#' @export
describeTemplate <- function(template, print = TRUE, var = NULL) {
 
  if (!file.exists(template)) {
    stop(paste("Template file", template, "does not exist."))
  }
  xml <- readLines(template)
 
  info <- .parseTemplateInfo(xml)
  checklist <- .structuralChecklist(xml)
 
  # every ${VAR} actually present, using the same regex the real
  # substitution step (write_clone_to_xml) uses to find required variables
  matches <- unlist(regmatches(xml, gregexpr("\\$\\{([^}]+)\\}", xml)))
  found_vars <- unique(sub("\\$\\{([^}]+)\\}", "\\1", matches))
 
  # Dowser-special placeholders aren't free-form kwargs -- drop them here,
  # they're already covered by the structural checklist above.
  special <- c("DATA", "CLONE", "DATE", "TRAIT", "TRAIT_NAME", "NODES", "MRCA",
               "OPERATORS", "ROOTFREQS", "EMP_EQFREQS", "ROOT_TYPE_PROBABILITIES",
               "NODES_TYPE_INIT")
  free_vars <- sort(setdiff(found_vars, special))
 
  get_field <- function(v, field, default) {
    if (!is.null(info$vars[[v]])) info$vars[[v]][[field]] else default
  }
 
  vars_df <- data.frame(
    name        = free_vars,
    documented  = vapply(free_vars, function(v) !is.null(info$vars[[v]]), logical(1)),
    tier        = vapply(free_vars, get_field, character(1), field = "tier", default = "undocumented"),
    description = vapply(free_vars, get_field, character(1), field = "description", default = ""),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
 
  struct_df <- data.frame(
    element = vapply(checklist, function(x) x$name, character(1)),
    present = vapply(checklist, function(x) x$present, logical(1)),
    effect  = vapply(checklist, function(x) if (x$present) x$if_present else x$if_absent, character(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  # append any author @struct notes for elements they chose to annotate
  for (nm in names(info$structs)) {
    idx <- which(vapply(struct_df$element, function(e) grepl(nm, e, fixed = TRUE), logical(1)))
    if (length(idx) > 0) {
      struct_df$effect[idx[1]] <- paste0(struct_df$effect[idx[1]], " Author note: ", info$structs[[nm]])
    } else {
      # note about something the checklist doesn't track by name (e.g. a
      # non-Dowser-special XML element the author wants to flag) -- append
      # as its own row rather than dropping it silently
      struct_df <- rbind(struct_df, data.frame(
        element = nm, present = NA, effect = paste("Author note:", info$structs[[nm]]),
        stringsAsFactors = FALSE
      ))
    }
  }
 
  if (!is.null(var)) {
    # Dowser-special placeholders map onto one row of the structural
    # checklist each (TRAIT and TRAIT_NAME share a row, since Dowser
    # treats them as one combined check). DATA and CLONE aren't in the
    # checklist at all -- they're unconditionally required, not optional
    # features, so there's no present/absent variability to report on.
    special_checklist_map <- c(
      DATE = "DATE placeholder",
      TRAIT = "TRAIT / TRAIT_NAME placeholder",
      TRAIT_NAME = "TRAIT / TRAIT_NAME placeholder",
      NODES = "NODES placeholder",
      MRCA = "MRCA placeholder",
      OPERATORS = "OPERATORS placeholder",
      ROOTFREQS = "ROOTFREQS placeholder",
      EMP_EQFREQS = "EMP_EQFREQS placeholder",
      ROOT_TYPE_PROBABILITIES = "ROOT_TYPE_PROBABILITIES placeholder",
      NODES_TYPE_INIT = "NODES_TYPE_INIT placeholder"
    )

    if (var %in% c("DATA", "CLONE")) {
      present <- any(grepl(paste0("\\$\\{", var, "\\}"), xml))
      effect <- if (present) {
        "Always required by Dowser. Present in this template."
      } else {
        "Always required by Dowser, but NOT present in this template -- this template cannot be used by Dowser at all."
      }
      var_result <- list(name = var, kind = "structural", present = present, effect = effect)
    } else if (var %in% names(special_checklist_map)) {
      row <- struct_df[struct_df$element == special_checklist_map[[var]], , drop = FALSE]
      var_result <- list(name = var, kind = "structural", present = row$present[1], effect = row$effect[1])
    } else if (var %in% vars_df$name) {
      row <- vars_df[vars_df$name == var, , drop = FALSE]
      var_result <- list(name = var, kind = "free", present = TRUE,
                          documented = row$documented[1], tier = row$tier[1],
                          description = row$description[1])
    } else {
      var_result <- list(name = var, kind = "not_found", present = FALSE)
    }

    if (print) {
      cat(var, " (in ", basename(template), "):\n", sep = "")
      if (var_result$kind == "not_found") {
        cat("  Not found -- no ${", var, "} placeholder anywhere in this template.\n", sep = "")
      } else if (var_result$kind == "structural") {
        status <- if (var_result$present) "[present]" else "[absent] "
        cat("  ", status, " (Dowser-recognized, not a free-form variable) - ", var_result$effect, "\n", sep = "")
      } else {
        tag  <- if (var_result$documented) var_result$tier else "UNDOCUMENTED"
        desc <- if (var_result$documented) var_result$description
                else "(no dowser-template-info entry for this variable)"
        cat("  [", tag, "]  ", desc, "\n", sep = "")
      }
    }

    return(invisible(var_result))
  }
 
  result <- list(description = info$description, vars = vars_df, structural = struct_df)
 
  if (print) {
    cat(basename(template), "\n\n")
 
    if (!is.null(info$description) && nchar(info$description) > 0) {
      cat(info$description, "\n\n")
    } else {
      cat("(no @description provided for this template)\n\n")
    }
 
    cat("Structural elements (built-in Dowser behavior):\n")
    for (i in seq_len(nrow(struct_df))) {
      status <- if (is.na(struct_df$present[i])) "[note]   "
                else if (struct_df$present[i]) "[present]" else "[absent] "
      cat(" ", status, struct_df$element[i], "-", struct_df$effect[i], "\n")
    }
    cat("\n")
    n_doc <- sum(vars_df$documented)
    n_var <- nrow(vars_df)
    cat("Variables (", n_var, " found, ", n_doc, " documented):\n", sep = "")
    for (i in seq_len(n_var)) {
      tag  <- if (vars_df$documented[i]) vars_df$tier[i] else "UNDOCUMENTED"
      desc <- if (vars_df$documented[i]) vars_df$description[i]
              else "(no dowser-template-info entry for this variable)"
      cat("  ", vars_df$name[i], " [", tag, "]  ", desc, "\n", sep = "")
    }
    if (n_var == 0) {
      cat("  (no free-form variables found in this template)\n")
    } else if (n_doc < n_var) {
      cat("\n", n_var - n_doc, " of ", n_var,
          " variables are undocumented in this template's dowser-template-info block.\n", sep = "")
    } else {
      cat("\nAll ", n_var, " variables are documented in this template's dowser-template-info block.\n", sep = "")
    }
  }
 
  invisible(result)
}

#' Strip a dowser-template-info comment block out of the final output XML
#'
#' Called at the end of write_clone_to_xml(), after the doc block has
#' already served its purpose (the missing-variable error hint). BEAST2
#' would just ignore an XML comment anyway, so this isn't required for
#' correctness -- it's here so the .xml files a person actually runs and
#' inspects don't carry documentation meant for template authors, and so a
#' stray literal ${VAR} written in comment prose (see the convention's
#' warning in TemplateInfo.R) can never linger in a real run's output.
#'
#' @param    xml   character vector of template lines
#' @return   character vector with the block (if any) removed entirely
#' @keywords internal
.stripTemplateInfo <- function(xml) {
  start_idx <- grep("^\\s*<!--\\s*dowser-template-info\\s*$", xml)
  if (length(start_idx) == 0) {
    return(xml) # nothing to strip
  }
  start_idx <- start_idx[1]

  end_idx <- grep("^\\s*-->\\s*$", xml)
  end_idx <- end_idx[end_idx > start_idx]
  if (length(end_idx) == 0) {
    return(xml) # malformed block; leave it alone rather than guess how much to remove
  }
  end_idx <- end_idx[1]

  xml[-(start_idx:end_idx)]
}






