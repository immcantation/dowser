## Functions for constructing clonal germline sequences
# Based closely on CreateGermlines.py
# Download IMGT GENE-DB databases
# 
# Loads all reference germlines from an Immcantation-formatted IMGT database.
# 
#TODO: make auto-download or internal IMGT database
# param dir      directory to contain database
# param organism  vector of species to download: human, mouse, rat, rabbit, rhesus_monkey
# param timeout  max time allowed for each download (temporarily sets global timeout option)
# return directory of downloaded IMGT BCR and TCR reference files
# details 
# Recoding of this script: 
# https://bitbucket.org/kleinstein/immcantation/src/master/scripts/fetch_imgtdb.sh
downloadIMGT <- function(dir="imgt", organism=c("human", "mouse", "rat", "rabbit", "rhesus_monkey"),
                         timeout=300){
  
  old_timeout <- getOption("timeout")
  options(timeout = max(timeout, getOption("timeout")))
  
  species_list <- organism
  OUTDIR <- dir
  REPERTOIRE <- "imgt"
  DATE <- format(Sys.time(), "%Y.%m.%d")
  
  # Associative array (for BASH v3) where keys are species folder names and values are query strings
  # rabbit:Oryctolagus_cuniculus, rat:Rattus_norvegicus, rhesus_monkey:Macaca_mulatta
  SPECIES_QUERY=c("human"="Homo+sapiens",
                  "mouse"="Mus",
                  "rat"="Rattus+norvegicus",
                  "rabbit"="Oryctolagus+cuniculus",
                  "rhesus_monkey"="Macaca+mulatta")
  # Associative array (for BASH v3) with species name replacements
  SPECIES_REPLACE=c("human"="Homo sapiens/Homo_sapiens",
                    "mouse"="Mus musculus/Mus_musculus",
                    "rat"="Rattus norvegicus/Rattus_norvegicus",
                    "rabbit"="Oryctolagus cuniculus/Oryctolagus_cuniculus",
                    "rhesus_monkey"="Macaca mulatta/Macaca_mulatta")
  
  for(SPECIES in species_list){
    KEY <- SPECIES
    VALUE <- SPECIES_QUERY[SPECIES]
    REPLACE_VALUE <- SPECIES_REPLACE[SPECIES]
    
    print(paste("Downloading", KEY,"repertoires into",OUTDIR))
    
    # Create directories
    dir.create(file.path(OUTDIR, KEY, "vdj"), recursive=TRUE, showWarnings=FALSE)
    dir.create(file.path(OUTDIR, KEY, "vdj_aa"), recursive=TRUE, showWarnings=FALSE)
    dir.create(file.path(OUTDIR, KEY, "leader_vexon"), recursive=TRUE, showWarnings=FALSE)
    dir.create(file.path(OUTDIR, KEY, "leader"), recursive=TRUE, showWarnings=FALSE)
    dir.create(file.path(OUTDIR, KEY, "constant"), recursive=TRUE, showWarnings=FALSE)
    
    # Download functions
    download_and_process <- function(url, outfile, replace) {
      tmpfile <- paste0(outfile, ".tmp")
      download.file(url, tmpfile, quiet=TRUE)
      lines <- readLines(tmpfile, warn=F)
      # fasta between last pre and /pre
      pre_start <- max(grep("<pre>", lines))
      pre_end <- max(grep("</pre>", lines))
      data <- lines[(pre_start + 1):(pre_end - 1)]
      
      # replace names with spaces with unerlines
      split <- strsplit(replace,split="/")[[1]]
      sub <- gsub(split[1], split[2], data)
      writeLines(sub, outfile)
      file.remove(tmpfile)
    }
    
    # Download VDJ regions
    cat("|---- Ig\n")
    for (CHAIN in c("IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", "IGLV", "IGLJ")) {
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=7.1+", CHAIN, "&species=", VALUE)
      FILE_NAME <- file.path(OUTDIR, KEY, "vdj", paste0(REPERTOIRE,"_", KEY, "_", CHAIN, ".fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE)
    }
    
    # Download leader and V exon for Ig
    for (CHAIN in c("IGHV", "IGKV", "IGLV")) {
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=8.1+", CHAIN, "&species=", VALUE,"&IMGTlabel=L-PART1+V-EXON")
      FILE_NAME <- file.path(OUTDIR, KEY, "leader_vexon", paste0(REPERTOIRE,"_lv_", KEY, "_", CHAIN, ".fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE)
    }
    
    # V amino acid for Ig
    for (CHAIN in c("IGHV", "IGKV", "IGLV")) {
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=7.3+",CHAIN,"&species=",VALUE)
      FILE_NAME <- file.path(OUTDIR, KEY, "vdj_aa", paste0(REPERTOIRE,"_aa_", KEY, "_", CHAIN, ".fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE)
    }
    
    # VDJ TCR
    cat("|---- TCR\n")
    for (CHAIN in c("TRAV", "TRAJ", "TRBV", "TRBD", "TRBJ", "TRDV", "TRDD", "TRDJ", "TRGV", "TRGJ")) {
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=7.1+",CHAIN,"&species=",VALUE)
      FILE_NAME <- file.path(OUTDIR, KEY, "vdj", paste0(REPERTOIRE,"_", KEY, "_", CHAIN, ".fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE)
    }
    
    # Download leader and V exon for Ig
    for (CHAIN in c("TRAV", "TRBV", "TRDV", "TRGV")) {
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=8.1+",CHAIN,"&species=",VALUE,"&IMGTlabel=L-PART1+V-EXON")
      FILE_NAME <- file.path(OUTDIR, KEY, "leader_vexon", paste0(REPERTOIRE,"_lv_", KEY, "_", CHAIN, ".fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE)
    }
    
    # V amino acid for TCR
    for (CHAIN in c("TRAV", "TRBV", "TRDV", "TRGV")) {
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=7.3+",CHAIN,"&species=",VALUE)
      FILE_NAME <- file.path(OUTDIR, KEY, "vdj_aa", paste0(REPERTOIRE,"_aa_", KEY, "_", CHAIN, ".fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE)
    }
    
    
    # Download leaders
    cat("|- Spliced leader regions\n")
    
    # Download leader and V exon for Ig
    cat("|---- Ig\n")
    for (CHAIN in c("IGH", "IGK", "IGL")) {
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=8.1+",CHAIN,"V&species=",VALUE,"&IMGTlabel=L-PART1+L-PART2")
      FILE_NAME <- file.path(OUTDIR, KEY, "leader", paste0(REPERTOIRE,"_", KEY, "_", CHAIN, "L.fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE)
    }
    
    # Download leader and V exon for TCR
    cat("|---- TCR\n")
    for (CHAIN in c("TRA", "TRB", "TRG", "TRD")) {
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=8.1+",CHAIN,"V&species=",VALUE,"&IMGTlabel=L-PART1+L-PART2")
      FILE_NAME <- file.path(OUTDIR, KEY, "leader", paste0(REPERTOIRE,"_", KEY, "_", CHAIN, "L.fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE)
    }
    
    # Download constant regions
    cat("|- Spliced constant regions\n")
    cat("|---- Ig\n")
    for (CHAIN in c("IGHC", "IGKC", "IGLC")) {
      # IMGT does not have artificially spliced IGKC / IGLC for multiple species
      if(CHAIN == "IGHC"){
        QUERY <- 14.1
      }else{
        QUERY <- 7.5
      }
      
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=",QUERY,"+",CHAIN,"&species=",VALUE)
      FILE_NAME <- file.path(OUTDIR, KEY, "constant", paste0(REPERTOIRE,"_", KEY, "_", CHAIN, ".fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE)
    }
    
    cat("|---- TCR\n")
    for (CHAIN in c("TRAC", "TRBC", "TRGC", "TRDC")) {
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=14.1+",CHAIN,"&species=",VALUE)
      FILE_NAME <- file.path(OUTDIR, KEY, "constant", paste0(REPERTOIRE,"_", KEY, "_", CHAIN, ".fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE)
    }
  }
  
  # Write download info
  INFO_FILE <- file.path(OUTDIR, "IMGT.yaml")
  info <- c(
    "source:  https://www.imgt.org/genedb",
    paste0("date:    ",DATE),
    "species:",
    paste0("    - ",paste0(species_list,":",SPECIES_QUERY[species_list]))
  )
  writeLines(info, con=INFO_FILE)
  
  options(timeout = old_timeout)
  
  return(dir)
}


#' \code{writeFasta} Write a fasta file of sequences given a 
#' named list of sequences
#' @param    seqs      named list of sequences (output from \code{readFasta})
#' @param    file      FASTA file for output
#'
#' @return   File of FASTA formatted sequences
#' @export
writeFasta <- function(seqs, file){
  if(!is.null(seqs)){
    out <- paste0(">", names(seqs), "\n", seqs)
    writeLines(out, con=file)
  }else{
    file.create(file)
  }
}

# Format IMGT database to Igblast database
# 
# Cleans IMGT database and creates blast dbs for IgBlast
# 
# param imgt  directory containing IMGT sequences ("dir" in downloadIMGT)
# param outdir   output directory for IgBlast databases, will contain /fasta and /database
# param igblast  IgBlast home directory location
# param organism   organism species to use (must be human, mouse, or rhesus_monkey)
# return outdir
# details 
# Recoding of this script: 
# https://bitbucket.org/kleinstein/immcantation/src/master/scripts/imgt2igblast.sh
imgtToIgblast <- function(imgt, outdir, igblast, 
                          organism=c("human", "mouse", "rhesus_monkey")){
  
  OUTDIR <- outdir
  GERMDIR <- imgt
  OUTFASTA <- file.path(OUTDIR, "fasta")
  OUTDATA <- file.path(OUTDIR, "database")
  
  exec = file.path(igblast, "bin", "makeblastdb")
  exec <- path.expand(exec)
  if(file.access(exec, mode=1) == -1) {
    stop("The file ", exec, " cannot be executed.")
  }
  
  dir.create(file.path(OUTDIR), recursive=TRUE, showWarnings=FALSE)
  dir.create(file.path(OUTFASTA), recursive=TRUE, showWarnings=FALSE)
  dir.create(file.path(OUTDATA), recursive=TRUE, showWarnings=FALSE)
  
  species_list <- dir(GERMDIR)
  species_list <- species_list[species_list != "IMGT.yaml"]
  species_list <- species_list[species_list %in% organism]
  
  
  vdj <- dir(file.path(GERMDIR, species_list, c("vdj")))
  vdj <- vdj[grepl(".fasta$", vdj)]
  vdj_aa <- dir(file.path(GERMDIR, species_list, c("vdj_aa")))
  vdj_aa <- vdj_aa[grepl(".fasta$", vdj_aa)]
  constant <- dir(file.path(GERMDIR, species_list, c("constant")))
  constant <- constant[grepl(".fasta$", constant)]
  
  # Create fasta files of each species, chain and segment combination
  for(SPECIES in species_list){
    #dir(file.path(GERMDIR, SPECIES, "vdj"))
    for(CHAIN in c("IG","TR")){
      # VDJ nucleotides
      for(SEGMENT in c("V","D","J")){
        #F=$(echo imgt_${SPECIES}_${CHAIN}_${SEGMENT}.fasta | tr '[:upper:]' '[:lower:]')
        #cat ${GERMDIR}/${SPECIES}/vdj/imgt_${SPECIES}_${CHAIN}?${SEGMENT}.fasta > ${TMPDIR}/${F}
        # requires updated readFasta
        fastas <- list.files(file.path(GERMDIR, SPECIES, "vdj"))
        fastas <- fastas[grepl(paste0("imgt_", SPECIES, "_", CHAIN, ".", SEGMENT,".*.fasta$"), fastas)]
        f <- unlist(lapply(fastas, function(x)readFasta(file.path(GERMDIR, SPECIES, "vdj", x))))
        cf <- cleanSeqs(f)
        writeFasta(cf, file.path(OUTFASTA, tolower(paste0("imgt_", SPECIES, "_", CHAIN, "_", SEGMENT, ".fasta"))))
      }
      
      # C nucleotides
      #F=$(echo imgt_${SPECIES}_${CHAIN}_c.fasta | tr '[:upper:]' '[:lower:]')
      #cat ${GERMDIR}/${SPECIES}/constant/imgt_${SPECIES}_${CHAIN}?C.fasta > ${TMPDIR}/${F}
      fastas <- list.files(file.path(GERMDIR, SPECIES, "constant"))
      fastas <- fastas[grepl(paste0("imgt_", SPECIES, "_", CHAIN, ".C.*.fasta$"), fastas)]
      f <- unlist(lapply(fastas, function(x)readFasta(file.path(GERMDIR, SPECIES, "constant", x))))
      cf <- cleanSeqs(f)
      writeFasta(cf, file.path(OUTFASTA, tolower(paste0("imgt_", SPECIES, "_", CHAIN, "_C", ".fasta"))))
      
      # V amino acids
      #F=$(echo imgt_aa_${SPECIES}_${CHAIN}_v.fasta | tr '[:upper:]' '[:lower:]')
      #cat ${GERMDIR}/${SPECIES}/vdj_aa/imgt_aa_${SPECIES}_${CHAIN}?V.fasta > ${TMPDIR}/${F}
      fastas <- list.files(file.path(GERMDIR, SPECIES, "vdj_aa"))
      fastas <- fastas[grepl(paste0("imgt_aa_", SPECIES, "_", CHAIN, ".V.*.fasta$"), fastas)]
      f <- unlist(lapply(fastas, function(x)readFasta(file.path(GERMDIR, SPECIES, "vdj_aa", x))))
      cf <- cleanSeqs(f)
      writeFasta(cf, file.path(OUTFASTA, tolower(paste0("imgt_aa_", SPECIES, "_", CHAIN, "_V", ".fasta"))))
    }
  }
  
  fastas <- list.files(OUTFASTA, full.names=FALSE)
  NT_FILES <- fastas[!grepl("imgt_aa_", fastas)]
  AA_FILES <- fastas[grepl("imgt_aa_", fastas)]
  
  for(file in NT_FILES){
    stem <- strsplit(file, split="\\.")[[1]][1]
    command <- paste("-parse_seqids -dbtype nucl -in",
                     file.path(OUTFASTA, file), "-out", file.path(OUTDATA, stem))
    params <- list(exec, command, stdout="", stderr="")
    status <- tryCatch(do.call(base::system2, params), error=function(e){
      return(e)
    }, warning=function(w){
      return(w)
    })
    if(status != 0){
      cat(paste(exec, command, "\n"))
      warning(status, paste(exec, command, "\n"))
    }
  }
  
  for(file in AA_FILES){
    stem <- strsplit(file, split="\\.")[[1]][1]
    command <- paste("-parse_seqids -dbtype prot -in",
                     file.path(OUTFASTA, file), "-out", file.path(OUTDATA, stem))
    params <- list(exec, command, stdout="", stderr="")
    status <- tryCatch(do.call(base::system2, params), error=function(e){
      return(e)
    }, warning=function(w){
      return(w)
    })
    if(status != 0){
      cat(paste(exec, command, "\n"))
      warning(status, paste(exec, command, "\n"))
    }
  }
  return(OUTDIR)
}

# remove IMGT gaps, uppercase sequences, and use second delimiter for name
cleanSeqs <- function(seqs, rm_gaps=TRUE){
  if(length(seqs) == 0){
    return(NULL)
  }
  if(rm_gaps){
    seqs <- toupper(gsub("\\.","",seqs))
  }else{
    seqs <- toupper(seqs)
  }
  names(seqs) <- sapply(strsplit(names(seqs), split="\\|"), function(x)x[2])
  # if repeated IDs, use the first one
  # might be consistent with cleanIMGT.py?
  cseqs <- list()
  for(x in 1:length(seqs)){
    if(is.null(cseqs[[names(seqs)[x]]])){
      cseqs[[names(seqs)[x]]] <- seqs[[x]]
    }
  }
  cseqs
}

# Runs IgBlast on a fasta file
#  
# param file      Fasta file of Ig or TR sequences
# param igblast   Location of IgBlast program directory, containing /bin/igblastn
# param refs     Reference directory of sequences (see downloadIMGT)
# param igdata  Internal data directory for IgBlast (defaults to \code{igblast}/internal_data)
# param organism    Organism (human, mouse, rhesus_monkey)
# param domain_system Currently only IMGT supported
# param outfile       Name of AIRR rearrangement file (must end in TSV)
# param nproc     Nummber of threads to use
# param db_prefix   File prefix for reference fastas
# param locus     Ig or TR
# param set_igdata  Set IGDATA environment variable?
# param return    Return data.frame of output?
# param verbose   Print extra info?
# return AIRR-rearrangement formatted data frame
# details 
# Runs IgBlast, similar to AssignGenes.py in Changeo.
# Must have IgBlast downloaded, precompiled binaries recommended:
# https://ncbi.github.io/igblast/
# 
# Note for M1/M2 Mac OS, may be necessary to install IgBlast .dmg and 
# manually set igdata. Otherwise most issues stem from the internal_data
# directory location.
# 
assignGenes <- function(
    file,
    igblast,
    refs,
    igdata = NULL,
    organism = "human",
    domain_system = "imgt",
    outfile = NULL,
    nproc = 1,
    db_prefix="imgt",
    locus = "Ig",
    set_igdata=TRUE,
    return=TRUE,
    verbose=TRUE){
  
  if(!organism %in% c("human", "mouse", "rhesus_monkey")){
    stop(paste("organism must be either human, mouse, rhesus_monkey"))
  }
  
  if(!locus %in% c("Ig", "TR")){
    stop(paste("locus must be either Ig or TR"))
  }
  
  exec <- file.path(igblast, "bin", "igblastn")
  exec <- path.expand(exec)
  if(file.access(exec, mode=1) == -1) {
    stop("The file ", exec, " cannot be executed.")
  }
  
  if(is.null(igdata)){
    print(paste0("igdata not specified, using ",igblast," for IGDATA"))
    igdata <- igblast
  }
  if(!is.null(igdata) && set_igdata){    
    cat(paste0("Setting IGDATA to ", igdata,"\n"))
    id <- file.path(igdata, "internal_data")
    opt <- file.path(igdata, "optional_file")
    if(!dir.exists(id)){
      stop(paste("Error setting IGDATA, directory:", id, "does not exist. See ?assignGenes for help."))
    }
    if(!dir.exists(opt)){
      stop(paste("Error setting IGDATA, directory:", opt, "does not exist. See ?assignGenes for help."))
    }
    Sys.setenv(IGDATA = igdata)
  }
  
  db_dir <- file.path(refs, "database")
  fasta_dir <- file.path(refs, "fasta")
  if(!dir.exists(db_dir)){
    stop(paste(db_dir, "does not exist."))
  }
  if(!dir.exists(fasta_dir)){
    stop(paste(db_dir, "does not exist."))
  }
  
  if(!grepl("\\.fasta$", file) && !grepl("\\.fa$", file)){
    stop(paste(file, "is not a fasta file (at least, it doesn't end in .fasta or .fa)"))
  }
  
  if(is.null(outfile)){
    filestem <- gsub("\\.fasta$|\\.fa$","", file)
  }else{
    if(!grepl("\\.tsv$", outfile)){
      stop(outfile, " must be a .tsv file")
    }
    filestem <- gsub("\\.tsv", "", outfile)
  }
  
  blastout <- paste0(filestem, ".tsv")
  #airrout <- paste0(filestem, ".tsv")
  
  db_v <- file.path(db_dir, paste(db_prefix, organism, tolower(locus), "v", sep="_"))
  db_d <- file.path(db_dir, paste(db_prefix, organism, tolower(locus), "d", sep="_"))
  db_j <- file.path(db_dir, paste(db_prefix, organism, tolower(locus), "j", sep="_"))
  aux <- file.path(igdata, "optional_file", paste0(organism,"_gl.aux"))
  
  
  stem <- strsplit(file, split="\\.")[[1]][1]
  command <- paste(
    paste0("-germline_db_V ", db_v),
    paste0("-germline_db_D ", db_d),
    paste0("-germline_db_J ", db_j),
    paste0("-auxiliary_data ", aux),
    paste0("-domain_system ", domain_system," -ig_seqtype ",locus," -organism ", organism),
    #paste0("-outfmt '7 std qseq sseq btop'"),
    paste0("-outfmt 19"),
    paste0("-query ", file),
    paste0("-num_threads ", nproc),
    paste0("-out ", blastout)
  )
  if(verbose){
    cat(paste(exec, command))
  }
  params <- list(exec, command, stdout="", stderr="")
  status <- tryCatch(do.call(base::system2, params), error=function(e){
    return(e)
  }, warning=function(w){
    return(w)
  })
  if(status != 0){
    cat(paste(exec, command, "\n"))
    stop("error running IgBlast")
  }
  
  if(return){
    return(airr::read_rearrangement(blastout))
  }
}

# Add IMGT gaps to IgBlast output
#  
# param db      AIRR-rearrangment formatted dataframe outputted from \link{assignGenes}
# param gapdb   Root directory of gapped fasta files (\code{dir} in \link{downloadIMGT})
# param organism  Organism (human, mouse, rhesus_monkey)
# param locus   Ig or TR
# param gapped_d  Include IMGT gaps in D regions? Only applicable to 
# return AIRR-rearrangement formatted data frame in which sequence_alignment and germline_alignment
# have been updated with IMGT gaps
# details 
# Similar functionality to MakeDb.py in Change-O. 

addGaps <- function(db, gapdb, organism="human", locus="Ig", gapped_d=FALSE){
  
  if(!organism %in% c("human", "mouse", "rhesus_monkey")){
    stop(paste("organism must be either human, mouse, rhesus_monkey"))
  }
  
  if(!locus %in% c("Ig", "TR")){
    stop(paste("locus must be either Ig or TR"))
  }
  
  gap_files <- list.files(file.path(gapdb, organism, "vdj"))
  if(length(gap_files) == 0){
    stop(paste("No files found in",file.path(gapdb, organism, "vdj")))
  }
  vgap_files <- gap_files[grepl(paste0("_",toupper(locus),".V\\.fasta$"), gap_files)]
  gaps <- unlist(lapply(vgap_files, function(x)readFasta(file.path(gapdb, organism, "vdj", x))))
  gaps <- cleanSeqs(gaps, rm_gaps=FALSE)
  
  # IGHD3-10*02 has a gap :-(
  # need to account for that
  dgap_files <- gap_files[grepl(paste0("_",toupper(locus),".D\\.fasta$"), gap_files)]
  dgaps <- unlist(lapply(dgap_files, function(x)readFasta(file.path(gapdb, organism, "vdj", x))))
  dgaps <- cleanSeqs(dgaps, rm_gaps=FALSE)
  
  
  results <- dplyr::tibble()
  for(i in 1:nrow(db)){
    row <- db[i,]
    
    insertion_sites <- 0
    
    if(is.na(row$np1_length) || is.na(row$np1)){
      row$np1_length <- 0
      row$np1 <- ""
    }
    if(is.na(row$np2_length) || is.na(row$np2)){
      row$np2_length <- 0
      row$np2 <- ""
    }
    if(is.na(row$d_sequence_alignment) || is.na(row$d_germline_alignment)){
      row$d_sequence_alignment <- ""
      row$d_germline_alignment <- ""
    }
    if(is.na(row$j_sequence_alignment) || is.na(row$j_germline_alignment)){
      row$j_sequence_alignment <- ""
      row$j_germline_alignment <- ""
    }
    if(is.na(row$v_sequence_alignment)){
      warning(row$sequence_id, " v alignment not found")
      next
    }
    
    if(paste0(row$v_sequence_alignment, row$np1, row$d_sequence_alignment, row$np2, row$j_sequence_alignment) != row$sequence_alignment){
      warning(row$sequence_id, " alignments don't add up")
      next
    }
    
    vcall <- strsplit(row$v_call, split=",")[[1]][1]
    dcall <- strsplit(row$d_call, split=",")[[1]][1]
    jcall <- strsplit(row$j_call, split=",")[[1]][1]
    
    # add IMGT gaps
    # remove germline-relative insertations
    qv <- strsplit(as.character(row$v_sequence_alignment),"")[[1]]
    gv <- strsplit(as.character(row$v_germline_alignment),"")[[1]]
    
    # remove insertions relative to germline
    v_del_pos <- gv == "-"
    gv <- gv[!v_del_pos]
    qv <- qv[!v_del_pos]
    
    qj <- strsplit(as.character(row$j_sequence_alignment),"")[[1]]
    gj <- strsplit(as.character(row$j_germline_alignment),"")[[1]]
    
    # remove insertions relative to germline
    j_del_pos <- gj == "-"
    gj <- gj[!j_del_pos]
    qj <- qj[!j_del_pos]
    
    insertion_sites <- sum(v_del_pos) + sum(j_del_pos)
    
    qd <- strsplit(as.character(row$d_sequence_alignment),"")[[1]]
    gd <- strsplit(as.character(row$d_germline_alignment),"")[[1]]
    
    # remove insertions relative to germline
    d_del_pos <- gd == "-"
    gd <- gd[!d_del_pos]
    qd <- qd[!d_del_pos]
    
    # v gene with IMGT gaps
    gappedv <- strsplit(gaps[[vcall]], split="")[[1]]
    
    # positions without gaps in IMGT
    no_gap_pos <- 1:length(gappedv)
    no_gap_pos <- no_gap_pos[gappedv != "."]
    
    # positions without gaps present in alignment
    v_germ_gap <- rep(".", length(gappedv))
    v_quer_gap <- rep(".", length(gappedv))
    
    # add non gap characters to non-gap sites
    gstart <- row$v_germline_start
    no_gap_pos_g <- no_gap_pos[1:row$v_germline_end] #gapped germline is always full length
    no_gap_pos <- no_gap_pos[gstart:(gstart + length(gv) - 1)]
    
    v_germ_gap[no_gap_pos_g] <- gappedv[no_gap_pos_g]
    v_quer_gap[no_gap_pos] <- qv
    
    if(length(v_germ_gap[no_gap_pos]) != length(gv)){
      stop("Error processing gapped V segment: ", row$sequence_id)
    }
    if(length(v_quer_gap[no_gap_pos]) != length(qv)){
      stop("Error processing gapped V segment: ", row$sequence_id)
    }
    
    # remove trailing gaps
    max_char <- max(which(v_germ_gap != "."))
    v_germ_gap <- v_germ_gap[1:max_char]
    v_quer_gap <- v_quer_gap[1:max_char]
    
    # correct germline coordinate positions
    vgaps_added <- sum(v_quer_gap == ".")
    vpos_del <- sum(v_del_pos)
    jpos_del <- sum(j_del_pos)
    
    row$v_germline_end <- max(no_gap_pos_g)
    # always start at 1 when adding IMGT gaps
    row$v_germline_start <- 1
    
    if(!is.na(dcall) && gapped_d){
      # d gene with IMGT gaps
      gappedd <- strsplit(dgaps[[dcall]], split="")[[1]]
      
      # positions without gaps in IMGT
      no_gap_pos <- 1:length(gappedd)
      no_gap_pos <- no_gap_pos[gappedd != "."]
      
      # positions without gaps present in alignment
      d_germ_gap <- rep(".", length(gappedd))
      d_quer_gap <- rep(".", length(gappedd))
      
      # add non gap characters to non-gap sites
      gdstart <- row$d_germline_start
      no_gap_pos_gd <- no_gap_pos[gdstart:row$d_germline_end] #non-gap positions in germline 
      no_gap_pos <- no_gap_pos[gdstart:row$d_germline_end] #non-gap positions in query
      
      d_germ_gap[no_gap_pos_gd] <- gappedd[no_gap_pos_gd]
      d_quer_gap[no_gap_pos] <- qd
      
      if(length(d_germ_gap[no_gap_pos]) != length(gd)){
        stop("Error processing gapped D segment: ", row$sequence_id)
      }
      if(length(d_quer_gap[no_gap_pos]) != length(qd)){
        stop("Error processing gapped D segment: ", row$sequence_id)
      }
      
      # remove trailing gaps
      max_char <- max(which(d_germ_gap != "."))
      d_germ_gap <- d_germ_gap[1:max_char]
      d_quer_gap <- d_quer_gap[1:max_char]
      
      # remove trailing gaps
      min_char <- min(which(d_germ_gap != "."))
      if(min_char > 1){
        d_germ_gap <- d_germ_gap[(min_char):length(d_germ_gap)]
        d_quer_gap <- d_quer_gap[(min_char):length(d_quer_gap)]
      }
      
      dgaps_added <- sum(d_quer_gap == ".")
      dpos_del <- sum(d_del_pos)
      row$d_germline_end <- max(no_gap_pos_gd) + dgaps_added
      
    }else{
      d_germ_gap <- strsplit(row$d_germline_alignment, split="")[[1]]
      d_quer_gap <- strsplit(row$d_sequence_alignment, split="")[[1]]
      dgaps_added <- 0
      dpos_del <- 0
    }
    
    vdj_quer_gapped <- paste0(
      paste(v_quer_gap,collapse=""),
      row$np1, 
      paste(d_quer_gap,collapse=""),
      row$np2,
      paste(qj,collapse=""))
    
    vdj_germ_gapped <- paste0(
      paste(v_germ_gap,collapse=""),
      paste0(rep("N",row$np1_length),collapse=""),
      paste(d_germ_gap,collapse=""),
      paste0(rep("N",row$np2_length),collapse=""),
      paste(gj,collapse=""))
    
    row$sequence_alignment <- vdj_quer_gapped
    row$germline_alignment <- vdj_germ_gapped
    
    # check region values
    # https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html
    if(!is.na(row$fwr1)){
      row$fwr1_start <- 1
      row$fwr1_end <- 26*3
      if(gsub("\\.|-","",substr(row$sequence_alignment, row$fwr1_start, row$fwr1_end)) != row$fwr1){
        if(insertion_sites == 0){
          warning(row$sequence_id, " fwr1 doesn't match")
          #next
        }else{
          warning(row$sequence_id, " fwr1 doesn't match but insertions present")
        }
      }#else{
      row$fwr1 <- substr(row$sequence_alignment, row$fwr1_start, row$fwr1_end)
      #}
    }
    if(!is.na(row$cdr1)){
      row$cdr1_start <- 26*3+1
      row$cdr1_end <- 38*3
      if(gsub("\\.|-","",substr(row$sequence_alignment, row$cdr1_start, row$cdr1_end)) != row$cdr1){
        if(insertion_sites == 0){
          warning(row$sequence_id, " cdr1 doesn't match")
          #next
        }else{
          warning(row$sequence_id, " cdr1 doesn't match but insertions present")
        }
      }#else{
      row$cdr1 <- substr(row$sequence_alignment, row$cdr1_start, row$cdr1_end)
      #}
    }
    if(!is.na(row$fwr2)){
      row$fwr2_start <- 38*3+1
      row$fwr2_end <- 55*3
      if(gsub("\\.|-","",substr(row$sequence_alignment, row$fwr2_start, row$fwr2_end)) != row$fwr2){
        if(insertion_sites == 0){
          warning(row$sequence_id, " fwr2 doesn't match")
          #next
        }else{
          warning(row$sequence_id, " fwr2 doesn't match but insertions present")
        }
      }#else{
      row$fwr2 <- substr(row$sequence_alignment, row$fwr2_start, row$fwr2_end)
      #}
    }
    if(!is.na(row$cdr2)){
      row$cdr2_start <- 55*3+1
      row$cdr2_end <- 65*3
      if(gsub("\\.|\\-","",substr(row$sequence_alignment, row$cdr2_start, row$cdr2_end)) != row$cdr2){
        if(insertion_sites == 0){
          warning(row$sequence_id, " cdr2 doesn't match")
          #next
        }else{
          warning(row$sequence_id, " cdr2 doesn't match but insertions present")
        }
      }#else{
      row$cdr2 <- substr(row$sequence_alignment, row$cdr2_start, row$cdr2_end)
      #}
    }
    if(!is.na(row$fwr3)){
      row$fwr3_start <- 65*3+1
      row$fwr3_end <- 104*3
      if(gsub("\\.|-","",substr(row$sequence_alignment, row$fwr3_start, row$fwr3_end)) != row$fwr3){
        if(insertion_sites == 0){
          warning(row$sequence_id, " fwr3 doesn't match")
          #next
        }else{
          warning(row$sequence_id, " fwr3 doesn't match but insertions present")
        }
      }#else{
      row$fwr3 <- substr(row$sequence_alignment, row$fwr3_start, row$fwr3_end)
      #}
    }
    if(!is.na(row$junction)){
      row$cdr3_start <- 104*3+1
      row$cdr3_end <- 104*3+1 + row$junction_length -6 - 1 + dgaps_added - dpos_del
      row$fwr4_start <- 104*3+1 + row$junction_length -6 - 1 + 1 + dgaps_added - dpos_del
      row$fwr4_end <- nchar(row$sequence_alignment)
      if(gsub("\\.","",substr(row$sequence_alignment, row$cdr3_start, row$cdr3_end)) != row$cdr3){
        if(insertion_sites == 0){
          warning(row$sequence_id, " cdr3 doesn't match")
          #next
        }else{
          warning(row$sequence_id, " cdr3 doesn't match but insertions present")
        }
      }#else{
      row$cdr3 <- gsub("\\.","",substr(row$sequence_alignment, row$cdr3_start, row$cdr3_end))
      #}
      #if(substr(row$sequence_alignment, row$fwr4_start, row$fwr4_end) != row$fwr4){
      # stop(row$sequence_id, " fwr4 doesn't match")
      #}
      row$fwr4 <- gsub("\\.","",substr(row$sequence_alignment, row$fwr4_start, row$fwr4_end))
      if(gsub("\\.","",substr(row$sequence_alignment, 310, 
                              310+row$junction_length-1 + dgaps_added - dpos_del)) != row$junction){
        if(insertion_sites == 0){
          warning(row$sequence_id, " junction doesn't match")
          next
        }else{
          warning(row$sequence_id, " junction doesn't - insertions present")
          next
        }
      }else{
        row$junction <- gsub("\\.","",substr(row$sequence_alignment, 310, 
                                             310+row$junction_length-1 + dgaps_added - dpos_del))
        row$junction_aa <- alakazam::translateDNA(row$junction)
      }
    }
    rm_cols <- c(
      "v_sequence_alignment", "d_sequence_alignment", "j_sequence_alignment",
      "v_germline_alignment", "d_germline_alignment", "j_germline_alignment",
      "v_alignment_start", "d_alignment_start", "j_alignment_start",
      "v_alignment_end", "d_alignment_end", "j_alignment_end",
      "v_sequence_alignment_aa", "d_sequence_alignment_aa", "j_sequence_alignment_aa",
      "v_germline_alignment_aa", "d_germline_alignment_aa", "j_germline_alignment_aa"
    )
    row <- row[,!names(row) %in% rm_cols]
    results <- dplyr::bind_rows(results, row)
  }
  return(results)
}

# Get the MRCA sequence of a particular tree
#
# param clone      A row of a \link{formatClones} object
# return The MRCA sequence of that clone
getMRCASeq <- function(clone){
  mrca_node <- ape::getMRCA(clone$trees[[1]], clone$trees[[1]]$tip.label)
  mrca <- clone$trees[[1]]$nodes[[mrca_node]]$sequence
  return(mrca)
}

# \code{updateAIRRGerm} Updates an AIRR-table to have a new germline alignment and associated germline statistics 
# 
# Uses the MRCA of a given clone as the germline alignment and updates the statistics by the associated IgBlast values
# 
# param airr_data      the airr-table used to create the clones object
# param clones   a clones object from \link{formatClones}
# param igblast  the file path to igblast
# param igblast_database the file path to the igblast database
# param references the file path to the references to be used with igblast 
# param organism the organism of the associated airr_data. Examples include "human" and "mouse"
# param locus if BCR data use "Ig" if TCR data use "TR"
# param outdir the directory to save intermediate files to
# param nproc the number of processors to use. Default is one. 
# return directory of downloaded IMGT BCR and TCR reference files
# 
updateAIRRGerm <- function(airr_data, clones, igblast, igblast_database, references, organism, locus, outdir, nproc = 1, ...){
  mrcas <- do.call(rbind, parallel::mclapply(1:nrow(clones), function(x){
    if(clones$data[[x]]@phylo_seq != "hlsequence"){
      value <- getMRCASeq(clones[x,])
      if(clones$data[[x]]@phylo_seq == "sequence"){
        seq_id <- paste0(clones$clone_id[x], "_heavy")
        og_germ <- airr_data$germline_alignment[airr_data$clone_id == clones$clone_id[x] &
                                                  airr_data$locus == "IGH"][1]
      } else{
        seq_id <- paste0(clones$clone_id[x], "_light")
        og_germ <- airr_data$germline_alignment[airr_data$clone_id == clones$clone_id[x] &
                                                  airr_data$locus != "IGH"][1]
      }
      numbers <- clones$data[[x]]@numbers
      gaps <- dplyr::setdiff(1:max(numbers), numbers)
      og_germ <- paste0(strsplit(og_germ, "")[[1]][-gaps], collapse = "")
      if(nchar(og_germ) > nchar(value)){
        diff <- nchar(og_germ) - nchar(value)
        pad_val <- substring(og_germ, nchar(og_germ) - diff + 1, nchar(og_germ))
        value <- paste0(value, pad_val)
      }
    }else{
      value <- getMRCASeq(clones[x,])
      restart_point <- which(diff(clones$data[[x]]@numbers) < 0) + 1
      hnumbers <- clones$data[[x]]@numbers[1:restart_point-1]
      lnumbers <- clones$data[[x]]@numbers[restart_point:length(clones$data[[x]]@numbers)]
      hvalue <- substring(value, 1, length(hnumbers))
      lvalue <- substring(value, length(hnumbers)+1, nchar(value))
      seq_id <- c(paste0(clones$clone_id[x], "_heavy"), paste0(clones$clone_id[x], "_light"))

      og_germ <- airr_data$germline_alignment[airr_data$clone_id == clones$clone_id[x] &
                                                airr_data$locus == "IGH"][1]
      gaps <- dplyr::setdiff(1:max(hnumbers), hnumbers)
      og_germ <- paste0(strsplit(og_germ, "")[[1]][-gaps], collapse = "")
      og_lgerm <- airr_data$germline_alignment[airr_data$clone_id == clones$clone_id[x] &
                                                 airr_data$locus != "IGH"][1]
      gaps <- dplyr::setdiff(1:max(lnumbers), lnumbers)
      og_lgerm <- paste0(strsplit(og_lgerm, "")[[1]][-gaps], collapse = "")
      if(nchar(og_germ) > nchar(hvalue)){
        diff <- nchar(og_germ) - nchar(hvalue)
        pad_val <- substring(og_germ, nchar(og_germ) - diff + 1, nchar(og_germ))
        hvalue <- paste0(hvalue, pad_val)
      }
      if(nchar(og_lgerm) > nchar(lvalue)){
        diff <- nchar(og_lgerm) - nchar(hvalue)
        pad_val <- substring(og_lgerm, nchar(og_lgerm) - diff + 1, nchar(og_lgerm))
        lvalue <- paste0(lvalue, pad_val)
      }
      value <- c(hvalue, lvalue)
    }
    temp <- data.frame(seq_id = seq_id,
                       sequence = value)
    return(temp)
  }, mc.cores = nproc))
  
  sink(path.expand(file.path(outdir, "mrca_seqs.fasta")))
  for(i in 1:nrow(mrcas)){
    cat(paste0(">", mrcas$seq_id[i]), "\n", mrcas$sequence[i], "\n", sep = "")
  }
  sink()
  
  ig_data <- assignGenes(file.path(outdir, "mrca_seqs.fasta"),
                         igblast = igblast,
                         refs = igblast_database,
                         outfile = file.path(outdir, "igblast.tsv"),
                         nproc = nproc, ...)
  ig_data <- addGaps(ig_data, gapdb=references, organism=organism, locus=locus, ...)
  ig_data$clone_id <- substring(ig_data$sequence_id, 1, nchar(ig_data$sequence_id) - 6)

  # get the length values 
  ig_data$v_germline_length <- ig_data$v_germline_end - ig_data$v_germline_start + 1
  ig_data$d_germline_length <- ig_data$d_germline_end - ig_data$d_germline_start + 1
  ig_data$j_germline_length <- ig_data$j_germline_end - ig_data$j_germline_start + 1
  
  # make the germline__mask column 
  ig_data$germline_alignment_d_mask <- unlist(parallel::mclapply(1:nrow(ig_data), function(x){
    vseq <- substring(ig_data$germline_alignment[x], 1, ig_data$v_germline_length[x])
    jseq <- substring(ig_data$germline_alignment[x], nchar(ig_data$germline_alignment[x]) - ig_data$j_germline_length[x] + 1, 
                      nchar(ig_data$germline_alignment[x]))
    nseq <- paste0(rep("N", nchar(ig_data$germline_alignment[x]) - sum(nchar(vseq), nchar(jseq))),
                    collapse = "")
    gmask <- paste0(vseq, nseq, jseq)
    return(gmask)
  }, mc.cores = nproc))
  write.table(ig_data, file = file.path(outdir, "igblast.tsv"), sep = "\t", row.names = FALSE)
  
  airr_data <- do.call(rbind, parallel::mclapply(1:nrow(airr_data), function(x){
    sub <- airr_data[x,]
    indx <- which(sub$clone_id == ig_data$clone_id & sub$locus == ig_data$locus)
    if(length(indx) > 0){
      sub$germline_alignment <- ig_data$germline_alignment[indx]
      sub$germline_alignment_d_mask <- ig_data$germline_alignment_d_mask[indx]
      sub$v_germline_start <- ig_data$v_germline_start[indx]
      sub$v_germline_end <- ig_data$v_germline_end[indx]
      sub$d_germline_start <- ig_data$d_germline_start[indx]
      sub$d_germline_end <- ig_data$d_germline_end[indx]
      sub$j_germline_start <- ig_data$j_germline_start[indx]
      sub$j_germline_end <- ig_data$j_germline_end[indx]
      sub$np1_length <- ig_data$np1_length[indx]
      sub$np2_length <- ig_data$np2_length[indx]
      sub$v_germline_length <- ig_data$v_germline_length[indx]
      sub$d_germline_length <- ig_data$d_germline_length[indx]
      sub$j_germline_length <- ig_data$j_germline_length[indx]
      sub$v_call <- ig_data$v_call[indx]
      sub$d_call <- ig_data$d_call[indx]
      sub$j_call <- ig_data$j_call[indx]
    }
    return(sub)
  }, mc.cores = nproc))
  
  # If you want to keep airr_data's original columns and update only the relevant ones:
  cols_to_update <- c("germline_alignment", "germline_alignment_d_mask", "v_germline_start", "v_germline_end",
                      "d_germline_start", "d_germline_end", "j_germline_start", "j_germline_end",
                      "np1_length", "np2_length", "v_germline_length", "d_germline_length", "j_germline_length",
                      "v_call", "d_call", "j_call")
  
  airr_data <- merge(airr_data, ig_data[, c("clone_id", "locus", cols_to_update)], 
                  by = c("clone_id", "locus"), all.x = TRUE, suffixes = c("", ".ig"))
  
  for (col in cols_to_update) {
    airr_data[[col]] <- airr_data[[paste0(col, ".ig")]]
    airr_data[[paste0(col, ".ig")]] <- NULL
  }
  return(airr_data)
}

#' \code{readIMGT} read in IMGT database
#' 
#' Loads all reference germlines from an Immcantation-formatted IMGT database.
#' 
# TODO: make auto-download or internal IMGT database
#' @param dir      directory containing Immcantation-formatted IMGT database
#' @param quiet    print warnings?
#' @return List of lists, leading to IMGT-gapped nucleotide sequences.
#' Structure of object is list[[locus]][[segment]]
#' locus refers to locus (e.g. IGH, IGK, TRA)
#' segment refers to gene segment category (V, D, or J)
#' @details Input directory must be formatted to Immcantation standard.
#' See https://changeo.readthedocs.io/en/stable/examples/igblast.html for example
#' of how to download.
#' @examples
#' # vdj_dir contains a minimal example of reference germlines 
#' # (IGHV3-11*05, IGHD3-10*01 and IGHJ5*02)
#' # which are the gene assignments for ExampleDb[1,]
#' vdj_dir <- system.file("extdata", "germlines", "imgt", "human", "vdj", package="dowser")
#' imgt <- readIMGT(vdj_dir)
#' @export
readIMGT <- function(dir, quiet=FALSE){
  sequences <- 0
  database <- list()
  files <- list.files(dir, full.names=TRUE)
  files <- files[grepl("\\.fasta$|\\.fa$", files)]
  if(length(files) == 0){
    stop("No fasta files found in directory")
  }
  for(file in files){
    fasta_list <- readFasta(file)
    fasta_list <- unlist(lapply(fasta_list,
                                function(x)toupper(paste0(x,collapse=""))))
    
    info <- strsplit(gsub("\\.fasta","",file), split="_")[[1]]
    length <- length(info)
    if(length < 3){
      stop(paste("Improperly formatted input file name:",file))
    }
    locus <- info[length]
    segment <- substr(info[length],4,4)
    locus <- substr(info[length],1,3)
    
    #less efficient, but deals with duplicate names like CreateGermlines
    #which uses the last allele available for a given duplicate
    #this happens in IMGT mouse database, which has the same genes from
    #multiple mouse strains
    fasta <- c()
    duplicates <- c()
    for(n in names(fasta_list)){
      name <- alakazam::getAllele(n, strip_d=FALSE)
      if(name %in% names(fasta)){
        duplicates <- c(duplicates,name)
      }else{
        fasta[name] <- fasta_list[[n]]
      }
    }
    if(length(duplicates) > 0 && !quiet){
      warning(paste("Segment IDs not unique in",
                    file,"\n",paste(duplicates,collapse=",")))
    }
    if(!locus %in% names(database)){
      database[[locus]] <- list()
    }
    database[[locus]][[segment]] <- fasta
    sequences <- sequences + length(fasta)
  }
  print(paste("Read in",sequences,"from",length(files),"fasta files"))
  
  database
}

#   TODO: this is not generalized for non-IMGT gapped sequences!
#   Note: This is separated into three functions in CreateGermlines.py
#' \link{getGermline} get germline segment from specified receptor and segment
#' @param receptor       row from AIRR-table containing sequence of interest
#' @param references     list of reference segments. Must be specific to  
#'                       locus and segment
#' @param segment        Gene segment to search. Must be V, D, or J.
#' @param field          Column name for segment gene call (e.g. v_call)
#' @param germ_start     Column name of index of segment start within germline 
#'                       segment (e.g. v_germline_start)
#' @param germ_end       Similar to germ_start, but specifies end of segment 
#'                       (e.g. v_germline_end) 
#' @param germ_length    Similar to germ_start, but specifies length of segment
#'                       (e.g. v_germline_end)
#' @param germ_aa_start  Column name of index of segment start within germline 
#'                       segment in AA (if amino_acid=TRUE, e.g. v_germline_start)
#' @param germ_aa_length Similar to germ_start, but specifies length of segment
#'                       in AA (if amino_acid=TRUE, e.g. v_germline_end)
#' @param amino_acid     Perform reconstruction on amino acid sequence (experimental)
#' @return String of germline sequence from specified segment aligned with the 
#' sequence in the seq column of \code{receptor}.
getGermline <- function(receptor, references, segment, field, 
                        germ_start, germ_end, germ_length, germ_aa_start,germ_aa_length, 
                        amino_acid=FALSE){
  # Extract allele call
  gene <- alakazam::getAllele(receptor[[field]], strip_d=FALSE)
  
  # Get germline start and length
  if(!amino_acid){
    pad_char <- 'N'
    start <- receptor[[germ_start]]
    len <- receptor[[germ_length]]
  }else{
    pad_char <- 'X'
    start <- receptor[[germ_aa_start]]
    len <- receptor[[germ_aa_length]]
  }
  if(is.na(start)){
    start <- 1
  }
  if(is.na(len)){
    len <- 0
  }
  
  # Build segment germline sequence
  if(segment == "V" || segment == "J"){
    if(is.na(gene)){
      germ_seq <- paste(rep(pad_char,len),collapse="")
    }else if(gene %in% names(references)){
      seq <- references[gene]
      seq_sub <- substr(seq, start, nchar(seq))
      pad <- len - nchar(seq_sub)
      if(pad < 0){
        pad <- 0
      }
      germ_seq <- paste0(substr(seq,start,start + len -1), 
                         paste(rep(pad_char, pad),collapse=""))
    }else{
      germ_seq <- NA
    }
  }else if(segment == "D"){
    if(is.na(gene)){
      germ_seq <- ""
    }else if(gene %in% names(references)){
      seq <- references[gene]
      germ_seq <- substr(seq, start, start + len -1)
    }else{
      germ_seq <- NA
    }
  }else{
    stop(paste("Segment",segment,"not found"))
  }
  
  if(is.na(germ_seq)){
    warning(paste("Allele",gene,
                  "is not in the provided germline database."))
  }
  return(germ_seq)
}

#  Assemble full length germline sequence
#  Arguments:
#    receptor (changeo.Receptor.Receptor): Receptor object
#    v_seq (str): V segment sequence as a string
#    d_seq (str): D segment sequence as a string
#    j_seq (str): J segment sequence as a string
#    amino_acid (bool): if True use X for N/P regions and amino acid positional fields,
#                       otherwise use N and nucleotide fields.
#  Returns:
#    str: full germline sequence
#' \link{stitchVDJ} combines germline gene segments to a single string
#' @param receptor       row from AIRR-table containing sequence of interest
#' @param v_seq          germline V segment sequence from \link{getGermline}
#' @param d_seq          germline D segment sequence from \link{getGermline}
#' @param j_seq          germline J segment sequence from \link{getGermline}
#' @param np1_length     Column name in receptor specifying np1 segment length 
#'                       (e.g. np1_length)
#' @param np2_length     Column name in receptor specifying np2 segment length 
#'                        (e.g. np1_length)
#' @param np1_aa_length  Column name in receptor specifying np1 segment length 
#'                        in AA (if amino_acid=TRUE, e.g. np1_length)
#' @param np2_aa_length Column name in receptor specifying np2 segment length 
#'                        in AA (if amino_acid=TRUE, e.g. np1_length)
#' @param amino_acid  Perform reconstruction on amino acid sequence (experimental)
#' @return Full length germline VDJ sequence aligned with aligned with the 
#' sequence in the \code{seq} column of \code{receptor}.
stitchVDJ <- function(receptor, v_seq, d_seq, j_seq, 
                      np1_length="np1_length", np2_length="np2_length",
                      np1_aa_length="np1_aa_length", np2_aa_length="np2_aa_length",
                      amino_acid=FALSE){
  # Get N/P lengths
  if(!amino_acid){
    np_char <- 'N'
    np1_len <- receptor[[np1_length]]
    np2_len <- receptor[[np2_length]]
  }else{
    np_char <- 'X'
    np1_len <- receptor[[np1_aa_length]]
    np2_len <- receptor[[np2_aa_length]]
  }
  if(is.na(np1_len)){
    np1_len <- 0
  }
  if(is.na(np2_len)){
    np2_len <- 0
  }
  # Assemble pieces starting with V segment
  sequence <- v_seq
  sequence <- paste0(sequence, paste(rep(np_char, np1_len),collapse=""))
  sequence <- paste0(sequence, d_seq)
  sequence <- paste0(sequence, paste(rep(np_char, np2_len),collapse=""))
  sequence <- paste0(sequence, j_seq)
  
  return(sequence)
}

#  Assemble full length region encoding
#  Arguments:
#    receptor (changeo.Receptor.Receptor): Receptor object
#    v_seq (str): V segment germline sequence as a string
#    d_seq (str): D segment germline sequence as a string
#    j_seq (str): J segment germline sequence as a string
#    amino_acid (bool): if True use amino acid positional fields, otherwise use nucleotide fields.
#  Returns:
#    str: string defining germline regions
#' \link{stitchRegions} Similar to \link{stitchVDJ} but with segment IDs 
#' instead of nucleotides
#' @param receptor      row from AIRR-table containing sequence of interest
#' @param v_seq         germline V segment sequence from \link{getGermline}
#' @param d_seq         germline D segment sequence from \link{getGermline}
#' @param j_seq         germline J segment sequence from \link{getGermline}
#' @param np1_length    Column name in receptor specifying np1 segment length 
#'                       (e.g. np1_length)
#' @param np2_length    Column name in receptor specifying np2 segment length 
#'                       (e.g. np1_length)
#' @param n1_length     Column name in receptor specifying n1 segment length 
#'                       (experimental)
#' @param n2_length     Column name in receptor specifying n2 segment length 
#'                       (experimental)
#' @param p3v_length    Column name in receptor specifying p3v segment length 
#'                       (experimental)
#' @param p5d_length    Column name in receptor specifying p5d segment length 
#'                       (experimental)
#' @param p3d_length    Column name in receptor specifying p3d segment length 
#'                       (experimental)
#' @param p5j_length    Column name in receptor specifying p5j segment length 
#'                       (experimental)
#' @param n2_length     Column name in receptor specifying n2 segment length 
#'                        (experimental)
#' @param np1_aa_length Column name in receptor specifying np1 segment length 
#'                        in AA (if amino_acid=TRUE, e.g. np1_length)
#' @param np2_aa_length Column name in receptor specifying np2 segment length 
#'                        in AA (if amino_acid=TRUE, e.g. np1_length)
#' @param amino_acid  Perform reconstruction on amino acid sequence (experimental)
#' @return Full length germline VDJ sequence with segment IDs instead of 
#' nucleotides.
#' @seealso \link{stitchVDJ}
stitchRegions <- function(receptor, v_seq, d_seq, j_seq, 
                          np1_length="np1_length", np2_length="np1_length",
                          n1_length="n1_length", p3v_length="p3v_length",
                          p5d_length="p5d_length", p3d_length="p3d_length",
                          n2_length="n2_length",p5j_length="p5j_length",
                          np1_aa_length="np1_aa_length", np2_aa_length="np2_aa_length",
                          amino_acid=FALSE){
  
  # Set mode for region definitions
  if(!is.null(receptor[[n1_length]])){
    full_junction <- TRUE
  }else{
    full_junction <- FALSE
  }
  
  # For now, don't support full_junction
  full_junction <- FALSE
  
  # Assemble pieces starting with V segment
  regions <- paste(rep('V',nchar(v_seq)),collapse="")
  
  # NP nucleotide additions after V
  if(amino_acid){
    # PNP nucleotide additions after V
    np1_len <- receptor[[np1_aa_length]]
    if(is.na(np1_len)){
      np1_len <- 0
    }
    regions <- paste0(regions, 
                      paste(rep('N', np1_len), collapse=""))
  }else if(!full_junction){
    # PNP nucleotide additions after V
    np1_len <- receptor[[np1_length]]
    if(is.na(np1_len)){
      np1_len <- 0
    }
    regions <- paste0(regions, 
                      paste(rep('N', np1_len), collapse=""))
  }else{
    # P nucleotide additions before N1
    p3v_len <- receptor[[p3v_length]]
    n1_len <- receptor[[n1_length]]
    p5d_len <- receptor[[p5d_length]]
    if(is.na(p3v_len)){
      p3v_len <- 0
    }
    if(is.na(n1_len)){
      n1_len <- 0
    }
    if(is.na(p5d_len)){
      p5d_len <- 0
    }
    
    # Update regions
    regions <- paste0(regions,paste(rep('P',p3v_len),collapse=""))
    regions <- paste0(regions,paste(rep('N',n1_len),collapse=""))
    regions <- paste0(regions,paste(rep('P',p5d_len),collapse=""))
  }
  # Add D segment
  regions <- paste0(regions, paste(rep('D', 
                                       nchar(d_seq)),collapse=""))
  
  # NP nucleotide additions before J
  if(amino_acid){
    np2_len <- receptor[[np2_aa_length]]
    if(is.na(np2_len)){
      np2_len <- 0
    }
    regions <- paste0(regions, 
                      paste(rep('N', np2_len), collapse=""))
  }else if(!full_junction){
    np2_len <- receptor[[np2_length]]
    if(is.na(np2_len)){
      np2_len <- 0
    }
    regions <- paste0(regions, 
                      paste(rep('N', np2_len), collapse=""))
  }else{
    p3d_len <- receptor[[p3d_length]]
    n2_len <- receptor[[n2_length]]
    p5j_len <- receptor[[p5j_length]]
    if(is.na(p3d_len)){
      p3d_len <- 0
    }
    if(is.na(n2_len)){
      n2_len <- 0
    }
    if(is.na(p5j_len)){
      p5j_len <- 0
    }
    
    # Update regions
    regions <- paste0(regions,paste(rep('P',p3d_len),collapse=""))
    regions <- paste0(regions,paste(rep('N',n2_len),collapse=""))
    regions <- paste0(regions,paste(rep('P',p5j_len),collapse=""))
  }
  # Add J segment
  regions <- paste0(regions, paste(rep('J', 
                                       nchar(j_seq)),collapse=""))
  
  return(regions)
}


#' \code{buildGermline} reconstruct germline segments from alignment data
#' 
#' Reconstruct germlines from alignment data.
#' 
#' @param receptor      row from AIRR-table containing sequence of interest
#' @param references    list of reference segments. Must be specific to locus
#' @param seq           Column name for sequence alignment
#' @param id            Column name for sequence ID
#' @param clone         Column name for clone ID
#' @param v_call        Column name for V gene segment gene call
#' @param d_call        Column name for D gene segment gene call
#' @param j_call        Column name for J gene segment gene call
#' @param v_germ_start  Column name of index of V segment start within germline 
#' @param v_germ_end    Column name of index of V segment end within germline 
#' @param v_germ_length Column name of index of V segment length within germline
#' @param d_germ_start  Column name of index of D segment start within germline 
#' @param d_germ_end    Column name of index of D segment end within germline 
#' @param d_germ_length Column name of index of D segment length within germline
#' @param j_germ_start  Column name of index of J segment start within germline 
#' @param j_germ_end    Column name of index of J segment end within germline 
#' @param j_germ_length Column name of index of J segment length within germline
#' @param np1_length    Column name in receptor specifying np1 segment length 
#' @param np2_length    Column name in receptor specifying np2 segment length
#' @param amino_acid    Perform reconstruction on amino acid sequence (experimental)
#' @return List of reconstructed germlines
#' @details Return object contains multiple IMGT-gapped germlines:
#' \itemize{
#'   \item  \code{full}:    Full length germline
#'   \item  \code{dmask}:   Full length germline with D region masked
#'   \item  \code{vonly}:   V gene segment of germline
#'   \item  \code{regions}: String showing VDJ segment of each position
#' }
#' @seealso \link{buildClonalGermline}, \link{stitchVDJ}
buildGermline <- function(receptor, references, 
                          seq="sequence_alignment", id="sequence_id", clone="clone_id",
                          v_call="v_call", d_call="d_call", j_call="j_call",
                          v_germ_start="v_germline_start",v_germ_end="v_germline_end",v_germ_length="v_germline_length",
                          d_germ_start="d_germline_start",d_germ_end="d_germline_end",d_germ_length="d_germline_length",
                          j_germ_start="j_germline_start",j_germ_end="j_germline_end",j_germ_length="j_germline_length",
                          np1_length="np1_length", np2_length="np2_length",
                          amino_acid=FALSE){
  
  # Build V segment germline sequence
  germ_vseq <- getGermline(receptor, references$V, segment="V",
                           field=v_call, germ_start=v_germ_start, germ_end=v_germ_end,
                           germ_length=v_germ_length, amino_acid=amino_acid)
  
  # Build D segment germline sequence
  germ_dseq <- getGermline(receptor, references$D, segment="D",
                           field=d_call, germ_start=d_germ_start,germ_end=d_germ_end,
                           germ_length=d_germ_length, amino_acid=amino_acid)
  
  # Build J segment germline sequence
  germ_jseq <- getGermline(receptor, references$J, segment="J",
                           field=j_call, germ_start=j_germ_start,germ_end=j_germ_end,
                           germ_length=j_germ_length, amino_acid=amino_acid)
  
  # Stitch complete germlines
  if(!is.na(germ_vseq) & !is.na(germ_dseq) & !is.na(germ_jseq)){
    germ_seq <- stitchVDJ(receptor, germ_vseq, germ_dseq, germ_jseq, 
                          np1_length=np1_length, np2_length=np2_length, amino_acid=amino_acid)
    regions <- stitchRegions(receptor, germ_vseq, germ_dseq, germ_jseq,
                             np1_length=np1_length, np2_length=np2_length, amino_acid=amino_acid)
    
    if(nchar(receptor[[seq]]) == 0){
      stop(paste("Sequence is missing from the sequence field",
                 receptor[[clone]]))
    }
    
    len_check <- nchar(germ_seq) - nchar(receptor[[seq]])
    if(len_check != 0){
      stop(paste("Germline sequence differs from input sequence by",
                 len_check,"in clone", receptor[[clone]], ", discarding"))
    }
    
    # Define return germlines object
    if(amino_acid){
      pad_char <- "X"
    }else{
      pad_char <- "N"
    }
    
    germ_dmask <- paste0(substr(germ_seq, 1, nchar(germ_vseq)),
                         paste(rep(pad_char,
                                   nchar(germ_seq) - nchar(germ_vseq) - nchar(germ_jseq)),
                               collapse=""))
    germ_dmask <- paste0(germ_dmask, substr(germ_seq, nchar(germ_dmask) + 1,
                                            nchar(germ_seq)))
    
    len_check <- nchar(germ_dmask) - nchar(receptor[[seq]])
    if(len_check != 0){
      stop(paste("Germline dmask sequence differs from input sequence by",
                 len_check,"in clone", receptor[[clone]], ", discarding"))
    }
  }else{
    germ_seq = NA
    germ_vseq = NA
    germ_dmask = NA
    regions= NA
  }
  
  germlines <- list()
  germlines$full <- germ_seq
  germlines$dmask <- germ_dmask
  germlines$vonly <- germ_vseq
  germlines$regions <- regions
  
  return(germlines)
}

#' \code{buildClonalGermline} Determine consensus clone sequence and create germline for clone
#' 
#' Determine consensus clone sequence and create germline for clone
#' 
#' @param receptors        AIRR-table containing sequences from one clone
#' @param references       Full list of reference segments, see \link{readIMGT}
#' @param chain            chain in \code{references} being analyzed
#' @param use_regions       Return string of VDJ regions? (optional)
#' @param vonly            Return germline of only v segment?
#' @param seq              Column name for sequence alignment
#' @param id               Column name for sequence ID
#' @param clone            Column name for clone ID
#' @param v_call           Column name for V gene segment gene call
#' @param j_call           Column name for J gene segment gene call
#' @param j_germ_length    Column name of J segment length within germline
#' @param j_germ_aa_length Column name of J segment amino acid length (if amino_acid=TRUE)
#' @param amino_acid       Perform reconstruction on amino acid sequence (experimental)
#' @param ...              Additional arguments passed to \link{buildGermline}
#' @return Tibble with reconstructed germlines
#' @details Return object adds/edits following columns:
#' \itemize{
#'   \item  \code{seq}:  Sequences potentially padded  same length as germline
#'   \item  \code{germline_alignment}: Full length germline
#'   \item  \code{germline_alignment_d_mask}: Full length, D region masked
#'   \item  \code{vonly}:   V gene segment of germline if vonly=TRUE
#'   \item  \code{regions}: String of VDJ segment in position if use_regions=TRUE
#' }
#' @seealso \link{createGermlines} \link{buildGermline}, \link{stitchVDJ}
buildClonalGermline <- function(receptors, references, 
                                chain="IGH", use_regions=FALSE, vonly=FALSE,
                                seq="sequence_alignment", id="sequence_id", clone="clone_id",
                                v_call="v_call", j_call="j_call", j_germ_length="j_germline_length",
                                j_germ_aa_length= "j_germline_aa_length",amino_acid=FALSE,...){
  
  if(amino_acid){
    stop("Amino acid mode not yet supported")
  }
  
  # Create dictionaries to count observed V/J calls
  v_dict <- c()
  j_dict <- c()
  
  # Amino acid settings
  if(amino_acid){
    pad_char <- 'X'
  }else{
    pad_char <- "N"
  }
  
  # Find longest sequence in clone, as well as V/J calls
  # note - always uses "first" for v/j calls
  v_dict <- unlist(lapply(receptors[[v_call]],function(x)
    alakazam::getAllele(x, strip_d=FALSE)))
  j_dict <- unlist(lapply(receptors[[j_call]],function(x)
    alakazam::getAllele(x, strip_d=FALSE)))
  seq_len <- unlist(lapply(receptors[[seq]],function(x)
    nchar(x)))
  
  # Consensus V and J having most observations
  vcounts <- table(v_dict)
  jcounts <- table(j_dict)
  v_cons <- names(vcounts)[vcounts == max(vcounts)]
  j_cons <- names(jcounts)[jcounts == max(jcounts)]
  max_len <- max(seq_len)
  
  # Consensus sequence(s) with consensus V/J calls and longest sequence
  cons_index <- v_dict %in% v_cons & j_dict %in% j_cons & seq_len == max_len
  
  # Consensus sequence(s) with consensus V/J calls but not the longest sequence
  if(sum(cons_index) == 0){
    cons_index <- v_dict == v_cons & j_dict == j_cons
  }
  
  # Return without germline if no sequence has both consensus V and J call
  if(sum(cons_index) == 0){
    warning(paste("Clone",unique(receptors[[clone]]),
                  "no sequence found with both consensus V and J calls."))
    germlines  <- list()
    germlines$full <- NA
    germlines$dmask <- NA
    germlines$regions <- NA
    germlines$vonly <- NA
  }else{
    # Select consensus Receptor, resolving ties by alphabetical ordering of sequence id.
    # CreateGermlines.py always sorts ids as characters
    cons_id <- sort(as.character(receptors[cons_index,][[id]]))[1]
    cons <- receptors[receptors[[id]] == cons_id,]
    
    # Pad end of consensus sequence with gaps to make it the max length
    gap_length <- max_len - nchar(cons[[seq]])
    if(gap_length > 0){
      if(amino_acid){
        cons[[j_germ_aa_length]] <- cons[[j_germ_aa_length]] + gap_length  
      }else{
        cons[[j_germ_length]] <- cons[[j_germ_length]] + gap_length
      }  
      cons[[seq]] <- paste0(cons[[seq]],
                            paste0(rep(pad_char,gap_length),collapse=""))
    }
    
    # Update lengths padded to longest sequence in clone
    receptors[[seq]] <- unlist(lapply(1:nrow(receptors),
                                      function(x){
                                        l = max_len - nchar(receptors[[seq]][x])
                                        if(amino_acid){
                                          receptors[[j_germ_aa_length]][x] = receptors[[j_germ_aa_length]][x] + l
                                        }else{
                                          receptors[[j_germ_length]][[x]] = receptors[[j_germ_length]][[x]] + l
                                        }
                                        paste0(receptors[[seq]][x],
                                               paste0(rep(pad_char,l),collapse=""))
                                      }))
    
    sub_db <- references[[chain]]
    
    if(length(sub_db) == 0){
      stop(paste("Reference database for",chain,"is empty"))
    }
    
    # Stitch consensus germline
    germlines <- tryCatch(buildGermline(cons, references=sub_db, seq=seq, 
                                        v_call=v_call, j_call=j_call, j_germ_length=j_germ_length,
                                        amino_acid=amino_acid,...),error=function(e)e)
    if("error" %in% class(germlines)){
      warning(paste("Clone",unique(receptors[[clone]]),
                    "germline reconstruction error.\n",
                    germlines))
      germlines  <- list()
      germlines$full <- NA
      germlines$dmask <- NA
      germlines$regions <- NA
      germlines$vonly <- NA
    }
  }
  
  receptors$germline_alignment <- germlines$full
  receptors$germline_alignment_d_mask <- germlines$dmask
  if(use_regions){
    receptors$regions <- germlines$regions
  }
  if(vonly){
    receptors$germline_alignment_vonly <- germlines$vonly
  }
  return(receptors)
}


#' \link{createGermlines} Determine consensus clone sequence and create germline for clone
#' @param data          AIRR-table containing sequences from one clone
#' @param references    Full list of reference segments, see \link{readIMGT}
#' @param locus         Name of the locus column in the input data
#' @param trim_lengths  Remove trailing Ns from \code{seq} column if length different from germline?
#' @param force_trim    Remove all characters from sequence if different from germline? (not recommended)
#' @param nproc         Number of cores to use
#' @param na.rm         Remove clones with failed germline reconstruction?
#' @param seq           Column name for sequence alignment
#' @param id            Column name for sequence ID
#' @param clone         Column name for clone ID
#' @param v_call        Column name for V gene segment gene call
#' @param d_call        Column name for D gene segment gene call
#' @param j_call        Column name for J gene segment gene call
#' @param v_germ_start  Column name of index of V segment start within germline
#' @param v_germ_end    Column name of index of V segment end within germline
#' @param v_germ_length Column name of index of V segment length within germline
#' @param d_germ_start  Column name of index of D segment start within germline
#' @param d_germ_end    Column name of index of D segment end within germline
#' @param d_germ_length Column name of index of D segment length within germline
#' @param j_germ_start  Column name of index of J segment start within germline
#' @param j_germ_end    Column name of index of J segment end within germline
#' @param j_germ_length Column name of index of J segment length within germline
#' @param np1_length    Column name in receptor specifying np1 segment length 
#' @param np2_length    Column name in receptor specifying np2 segment length
#' @param amino_acid    Perform reconstruction on amino acid sequence (experimental)
#' @param fields        Character vector of additional columns to use for grouping. 
#'                      Sequences with disjoint values in the specified fields 
#'                      will be considered as separate clones.
#' @param verbose       amount of rubbish to print
#' @param ...           Additional arguments passed to \link{buildGermline}
#' @return Tibble with reconstructed germlines
#' @details Return object adds/edits following columns:
#' \itemize{
#'   \item  \code{seq}:  Sequences potentially padded  same length as germline
#'   \item  \code{germline_alignment}: Full length germline
#'   \item  \code{germline_alignment_d_mask}: Full length, D region masked
#'   \item  \code{vonly}:   V gene segment of germline if vonly=TRUE
#'   \item  \code{regions}: String of VDJ segment in position if use_regions=TRUE
#' }
#' @seealso \link{createGermlines} \link{buildGermline}, \link{stitchVDJ}
#' @examples 
#' vdj_dir <- system.file("extdata", "germlines", "imgt", "human", "vdj", package="dowser")
#' imgt <- readIMGT(vdj_dir)
#' db <- createGermlines(ExampleAirr[1,], imgt)
#' @export
createGermlines <- function(data, references, locus="locus", trim_lengths=FALSE, force_trim=FALSE,
                            nproc=1, seq="sequence_alignment", v_call="v_call", d_call="d_call", 
                            j_call="j_call", amino_acid=FALSE,  id="sequence_id", clone="clone_id",
                            v_germ_start="v_germline_start", v_germ_end="v_germline_end", v_germ_length="v_germline_length",
                            d_germ_start="d_germline_start", d_germ_end="d_germline_end", d_germ_length="d_germline_length",
                            j_germ_start="j_germline_start", j_germ_end="j_germline_end", j_germ_length="j_germline_length",
                            np1_length="np1_length", np2_length="np2_length", na.rm=TRUE, fields=NULL, verbose=0, ...){
  
  if(nrow(data) == 0){
    warning("No data provided!")
    return(data)
  }
  if(locus %in% c("IGH", "IGK", "IGL")){
    stop(paste0("locus option now indicates locus column name, not value. Sorry for the change!",
                " createGermlines now does all loci at once, so no need to separate by locus."))
  }
  if(!locus %in% names(data)){
    warning(paste0(locus, " column not found, attempting to extract locus from V call"))
    data[[locus]] = substr(data[[v_call]],1,3)
    warning(paste("Loci found:",unique(data[[locus]])))
  }
  complete <- dplyr::tibble()
  required <- c(seq, id, clone, 
                np1_length, np1_length, 
                v_call, d_call, j_call,
                v_germ_start, v_germ_end,
                d_germ_start, d_germ_end,
                j_germ_start, j_germ_end, locus, fields)
  if(sum(!required %in% names(data)) != 0){
    stop(paste("Required columns not found in data:",
               paste(required[!required %in% names(data)],collapse=", ")))
  }
  if(sum(is.na(data[[clone]])) > 0){
    stop("NA values in clone id column found, please remove.")
  }
  
  # check if there are "" in the d_call column instead of NAs CGJ 11/1/23
  data[[d_call]][data[[d_call]] == ""] <- NA
  has_dup_ids <- max(table(data %>% select(!!!rlang::syms(c(id, fields))))) != 1
  if (has_dup_ids){
    stop("Sequence IDs are not unique!")
  }

  if(!v_germ_length %in% names(data)){
    data[[v_germ_length]] <- data[[v_germ_end]] - data[[v_germ_start]] + 1
  }
  if(!d_germ_length %in% names(data)){
    data[[d_germ_length]] <- data[[d_germ_end]] - data[[d_germ_start]] + 1
  }
  if(!j_germ_length %in% names(data)){
    data[[j_germ_length]] <- data[[j_germ_end]] - data[[j_germ_start]] + 1
  }
  if(sum(is.na(data[[v_germ_length]])) > 0){
    data[[v_germ_length]][is.na(data[[v_germ_length]])] = 
      data[[v_germ_end]][is.na(data[[v_germ_length]])] -
      data[[v_germ_start]][is.na(data[[v_germ_length]])] + 1
  }
  if(sum(is.na(data[[d_germ_length]])) > 0){
    data[[d_germ_length]][is.na(data[[d_germ_length]])] = 
      data[[d_germ_end]][is.na(data[[d_germ_length]])] -
      data[[d_germ_start]][is.na(data[[d_germ_length]])] + 1
  }
  if(sum(is.na(data[[j_germ_length]])) > 0){
    data[[j_germ_length]][is.na(data[[j_germ_length]])] = 
      data[[j_germ_end]][is.na(data[[j_germ_length]])] -
      data[[j_germ_start]][is.na(data[[j_germ_length]])] + 1
  }
  
  if(sum(is.na(data[[v_germ_length]])) > 0 | 
     sum(is.na(data[[j_germ_length]])) > 0){
    stop("Missing values in v_germ_length or j_germ_length")
  }


  # check if sequence_alignments contain trailing Ns and trim if desired
  # trailing Ns frequently cause length errors downstream
  # KBH 8/5/24
  g_lengths <- sapply(1:nrow(data), function(x)sum(data[[v_germ_length]][x], data[[np1_length]][x], 
    data[[d_germ_length]][x], data[[np2_length]][x], data[[j_germ_length]][x], na.rm=TRUE))
  g_diffs <- nchar(data[[seq]]) - g_lengths
  if(sum(g_diffs > 0) > 0){
    if(!trim_lengths && !force_trim){
     warning(sum(g_diffs)," sequence lengths longer than predicted germlines, consider setting ",
      "trim_lengths=TRUE if germlines fail")
    }else{
      too_short <- data[g_diffs > 0,]
      too_short_diffs <- g_diffs[g_diffs > 0]
      too_short_starts <- nchar(too_short[[seq]]) - too_short_diffs
      # short_seqs <- strsplit(too_short[[seq]], split="")
      to_cut <- sapply(1:nrow(too_short), function(x){
        substr(too_short[[seq]][x],too_short_starts[x] + 1, nchar(too_short[[seq]][x]))
      })
      atcg <- grepl("[ATCG]",to_cut)
      too_short[[seq]][!atcg] <- sapply(1:nrow(too_short[!atcg,]), function(x){
        substr(too_short[[seq]][!atcg][x], 1, too_short_starts[!atcg][x])
      })
      cat("Trimmed ",sum(!atcg),
        "sequences that differed from predicted germline only by non-ATCG characters.",
        sum(atcg), "differed by ATCG characters.\n")
      if(sum(atcg) > 0 && !force_trim){
        cat("Can remove ATCG characters if force_trim=TRUE, but this may indicate misalignment of you data.\n")
      }
      if(force_trim){
        cat("Forcibly removing ATCG characters from", sum(atcg), "sequences\n")
        too_short[[seq]][atcg] <- sapply(1:nrow(too_short[atcg,]), function(x){
          substr(too_short[[seq]][atcg][x], 1, too_short_starts[atcg][x])
        })
      }
      m <- match(data[[id]], too_short[[id]])
      seqs <- too_short[[seq]][m]
      seqs[is.na(m)] <- data[[seq]][is.na(m)]
      data[[seq]] <- seqs
    }
  }
  unique_clones <- unique(data[,unique(c(clone,fields)),drop=F])
  data[['tmp_row_id']] <- 1:nrow(data)
  complete <- parallel::mclapply(1:nrow(unique_clones), function(x){
    sub <- dplyr::right_join(data, unique_clones[x,,drop=F], by=c(clone,fields))
    if(verbose > 0){
      print(unique(sub[[clone]]))
    }
    glines <- lapply(unique(sub[[locus]]), function(l){
      buildClonalGermline(
        sub[sub[[locus]] == l,], 
        references=references,
        chain=l,
        seq=seq,
        v_call=v_call,
        d_call=d_call,
        j_call=j_call,
        amino_acid=amino_acid,
        id =id ,
        clone=clone,
        v_germ_start=v_germ_start,
        v_germ_end=v_germ_end,
        v_germ_length=v_germ_length,
        d_germ_start=d_germ_start,
        d_germ_end=d_germ_end,
        d_germ_length=d_germ_length,
        j_germ_start=j_germ_start,
        j_germ_end=j_germ_end,
        j_germ_length=j_germ_length,
        np1_length=np1_length,
        np2_length=np2_length,
        ...)
    })
    gline <- dplyr::bind_rows(glines)
    gline
  }, mc.cores=nproc)
  results <- dplyr::bind_rows(complete) %>%
    arrange(!!rlang::sym("tmp_row_id")) %>%
    select(-!!rlang::sym("tmp_row_id"))
  if(na.rm){
    bad_clones <- unique(results[is.na(results$germline_alignment_d_mask),][[clone]])
    if(dplyr::n_distinct(bad_clones) > 0){
      warning(paste("Removing",
                    dplyr::n_distinct(bad_clones),"failed clonal germlines. Clones:",
                    paste(bad_clones,collapse=",")))
      results <- results[!is.na(results$germline_alignment_d_mask),]
    }
  }
  results
}

# finds what the consensus sequence was that createGermlines used.
findConsensus <- function(receptors, clone_id, v_call = "v_call", j_call = "j_call",
                          seq = "sequence_alignment", id = "sequence_id"){
  v_dict <- c()
  j_dict <- c()
  pad_char <- "N"
  v_dict <- unlist(lapply(receptors[[v_call]],function(x)
    alakazam::getAllele(x, strip_d=FALSE)))
  j_dict <- unlist(lapply(receptors[[j_call]],function(x)
    alakazam::getAllele(x, strip_d=FALSE)))
  seq_len <- unlist(lapply(receptors[[seq]],function(x)
    nchar(x)))
  
  # Consensus V and J having most observations
  vcounts <- table(v_dict)
  jcounts <- table(j_dict)
  v_cons <- names(vcounts)[vcounts == max(vcounts)]
  j_cons <- names(jcounts)[jcounts == max(jcounts)]
  max_len <- max(seq_len)
  
  # Consensus sequence(s) with consensus V/J calls and longest sequence
  cons_index <- v_dict %in% v_cons & j_dict %in% j_cons & seq_len == max_len
  
  # Consensus sequence(s) with consensus V/J calls but not the longest sequence
  if(sum(cons_index) == 0){
    cons_index <- v_dict == v_cons & j_dict == j_cons
  }
  cons_id <- sort(as.character(receptors[cons_index,][[id]]))[1]
  cons <- receptors[receptors[[id]] == cons_id,]
  rec_v <- unlist(lapply(cons[[v_call]],function(x)
    alakazam::getAllele(x, strip_d=FALSE)))
  rec_j <- unlist(lapply(cons[[j_call]],function(x)
    alakazam::getAllele(x, strip_d=FALSE)))
  temp <- data.frame(clone_id = clone_id, locus = receptors$locus[1], v_call = rec_v, 
                     j_call = rec_j, cons_id = cons_id)
  return(temp)
}

# checkGenesUCA is what is run if check_genes = TRUE in getTreesAndUCAs
# sub is an clones object for only 1 clone
# data is the airr data used to formatclones
# v is a string for the v gene
# mracdr3 is a string of the junction (conserved site to conserved site)
# j is a string for the j gene 
# references is the readIMGT output 
# tree_df is the the codon table 
# subDir is where things should be saved
# chain is used to indicate heavy chain ("H") or light chain ("L")
checkGenesUCA <- function(sub, data, v, mrcacdr3, j, references, tree_df, subDir,
                          clone_ids, r, chain = "H"){
  if(chain == "H"){
    cons <- findConsensus(receptors = dplyr::filter(data, !!rlang::sym("clone_id") == sub$clone_id &
                                          !!rlang::sym('locus') == "IGH"), clone_id = sub$clone_id)
  } else{
    cons <- findConsensus(receptors = dplyr::filter(data, !!rlang::sym("clone_id") == sub$clone_id &
                                          !!rlang::sym('locus') != "IGH"), clone_id = sub$clone_id)
  }
  cons <- data[data$sequence_id == cons$cons_id,]
  if(sub$data[[1]]@phylo_seq == "hlsequence"){
    numbers <- sub$data[[1]]@numbers
    restart_point <- which(diff(numbers) < 0) + 1
    if(chain == "H"){
      numbers <- numbers[1:restart_point - 1]
    } else{
      numbers <- numbers[restart_point: length(numbers)]
    }
  }else{
    numbers <- sub$data[[1]]@numbers
  }
  gaps <- dplyr::setdiff(1:max(numbers), numbers)
  uca <- strsplit(paste0(v, mrcacdr3, j), "")[[1]]
  for(pos in gaps) {
    uca <- append(uca, ".", after = pos - 1)
  }
  uca <- paste(uca, collapse = "")
  for(pos in gaps) {
    r <- append(r, "gap", after = pos - 1)
  }
  ref_v <- references[[cons$locus]]$V[which(names(references[[cons$locus]]$V) == 
                                              strsplit(cons$v_call, ",")[[1]][1])]
  ref_v <- substring(ref_v, cons$v_germline_start, cons$v_germline_end)
  if(chain == "L" | sub$data[[1]]@phylo_seq == "lsequence"){
    pad_length <- sum(strsplit(substring(sub$data[[1]]@lgermline, (nchar(sub$data[[1]]@lgermline)-2),
                                         (nchar(sub$data[[1]]@lgermline))), "")[[1]] == "N")
  } else{
    pad_length <- sum(strsplit(substring(sub$data[[1]]@germline, (nchar(sub$data[[1]]@germline)-2),
                                         (nchar(sub$data[[1]]@germline))), "")[[1]] == "N")
  }
  ref_j <- references[[cons$locus]]$J[which(names(references[[cons$locus]]$J) == 
                                              strsplit(cons$j_call, ",")[[1]][1])]
  # do a check if nchar(germline alignment) is > sum(igblast stats) 
  # if greater add to j length 
  if(is.na(cons$d_germline_length)){
    cons$d_germline_length <- 0 
  } 
  if(is.na(cons$np2_length)){
    cons$np2_length <- 0 
  }
  igblast_len <- sum(as.numeric(cons$v_germline_length), as.numeric(cons$np1_length),
                     as.numeric(cons$d_germline_length), as.numeric(cons$np2_length),
                     as.numeric(cons$j_germline_length))
  if(igblast_len < nchar(cons$germline_alignment)){
    ig_diff <- nchar(cons$germline_alignment) - igblast_len
    cons$j_germline_length <- cons$j_germline_length + ig_diff
    cons$j_germline_end <- cons$j_germline_end + ig_diff
  }
  ref_j <- substring(ref_j, cons$j_germline_start, cons$j_germline_end)
  if(nchar(cons$germline_alignment) > nchar(uca)){
    # trim the j gene 
    diff <- nchar(cons$germline_alignment) - nchar(uca)
    ref_j <- substring(ref_j, 1, nchar(ref_j)-diff)
  }
  ref_j <- paste0(ref_j, paste(rep("N", pad_length), collapse = ""))
  ref_v <- paste0(strsplit(ref_v, "")[[1]][-gaps], collapse = "")
  j_gaps <- gaps[gaps >= nchar(uca) - nchar(ref_j) + 1]
  
  if(length(j_gaps) > 0){
    j_gaps <- j_gaps - (nchar(uca) - nchar(ref_j) + 1)
    ref_j <- paste0(strsplit(ref_j, "")[[1]][-j_gaps], collapse = "")
  }
  if(cons$locus == "IGH"){
    if(nchar(ref_j) > (nchar(sub$data[[1]]@germline) - sum(nchar(ref_v),  
                                                           cons$np1_length, cons$d_germline_length, cons$np2_length))){
      diff <- nchar(ref_j) - (nchar(sub$data[[1]]@germline) - sum(nchar(ref_v),  
                                                                  cons$np1_length, cons$d_germline_length, cons$np2_length))
      ref_j <- substring(ref_j, 1, nchar(ref_j) - diff)
    }
  } else{
    if(nchar(ref_j) > (nchar(sub$data[[1]]@lgermline) - sum(nchar(ref_v),  
                                                           cons$np1_length, cons$d_germline_length, cons$np2_length))){
      diff <- nchar(ref_j) - (nchar(sub$data[[1]]@lgermline) - sum(nchar(ref_v),  
                                                                  cons$np1_length, cons$d_germline_length, cons$np2_length))
      ref_j <- substring(ref_j, 1, nchar(ref_j) - diff)
    }
  }
  v_groups <- lapply(seq(1, nchar(ref_v), by = 3),
                     function(i) i:min(i+2, nchar(ref_v)))
  ref_v <- strsplit(ref_v, "")[[1]]
  if(nchar(ref_j) %% 3 == 0) {
    # Perfect multiple of 3, normal grouping
    j_groups <- lapply(seq(1, nchar(ref_j), by = 3),
                       function(i) i:min(i+2, nchar(ref_j)))
  } else{
    # Start with incomplete group, then complete groups of 3
    j_groups <- list()
    j_groups[[1]] <- 1:(nchar(ref_j) %% 3)  # First incomplete group
    
    # Add complete groups of 3
    start_positions <- seq((nchar(ref_j) %% 3) + 1, nchar(ref_j), by = 3)
    complete_groups <- lapply(start_positions,
                              function(i) i:min(i+2, nchar(ref_j)))
    
    j_groups <- c(j_groups, complete_groups)
  }
  ref_j <- strsplit(ref_j, "")[[1]]
  # update the tree_df to only have the values at each site for the v gene 
  v_cons <- (nchar(v) + 1): (nchar(v) + 3)
  v_con_indx <- which(sapply(v_groups, function(x) any(x %in% v_cons)))
  v_df <- do.call(rbind, lapply(1:length(v_groups), function(i){
    temp <- tree_df[tree_df$site == i - 1,]
    if(length(v_con_indx) > 0){
      if(i != v_con_indx){
        if(length(v_groups[[i]]) == 3){
          temp <- temp[temp$codon == paste0(ref_v[v_groups[[i]]], collapse = ""),]
        } else{
          values <- paste0(ref_v[v_groups[[i]]], collapse = "")
          temp <- temp[startsWith(temp$codon, values), ]
        }
      } else{
        if(length(v_groups[[i]]) == 3){
          temp <- temp[
            (temp$codon == paste0(ref_v[v_groups[[i]]], collapse = "")) |
              (alakazam::translateDNA(temp$codon) == "C"),]
          if(alakazam::translateDNA(paste0(ref_v[v_groups[[i]]], collapse = "")) == "C"){
            temp <- temp[temp$codon == paste0(ref_v[v_groups[[i]]], collapse = ""),]
          }
        } else{
          values <- paste0(ref_v[v_groups[[i]]], collapse = "")
          temp <- temp[(startsWith(temp$codon, values)) |
                         (alakazam::translateDNA(temp$codon) == "C"), ]
        }
      }
    } else{
      if(length(v_groups[[i]]) == 3){
        temp <- temp[temp$codon == paste0(ref_v[v_groups[[i]]], collapse = ""),]
      } else{
        values <- paste0(ref_v[v_groups[[i]]], collapse = "")
        temp <- temp[startsWith(temp$codon, values), ]
      }
    }
    return(temp)
  }))
  j_df <- tree_df[tree_df$site %in% tail(sort(unique(tree_df$site)),
                                         length(j_groups)), ]
  j_df$new_site <- j_df$site - (min(j_df$site) - 1)
  j_con <- (nchar(v)+ nchar(mrcacdr3) - 2): (nchar(v)+ nchar(mrcacdr3))
  if(chain == "L" | sub$data[[1]]@phylo_seq == "lsequence"){
    offset <- nchar(sub$data[[1]]@lgermline) - max(j_groups[[length(j_groups)]])
  } else{
    offset <- nchar(sub$data[[1]]@germline) - max(j_groups[[length(j_groups)]])
  }
  j_groups_num <- lapply(j_groups, function(x) x + offset)
  j_con_indx <- which(sapply(j_groups_num, function(x) any(x %in% j_con)))
  j_df <- do.call(rbind, lapply(1:length(j_groups), function(i){
    temp <- j_df[j_df$new_site == i,]
    if(length(j_con_indx) > 0){
      if(i != j_con_indx){
        if(length(j_groups[[i]]) == 3){
          if(sum("N" %in% ref_j[j_groups[[i]]]) == 0){
            if(alakazam::translateDNA(paste0(ref_j[j_groups[[i]]], collapse = "")) != "*"){
              temp <- temp[temp$codon == paste0(ref_j[j_groups[[i]]], collapse = ""),]
            }
          } else{
            values <- ref_j[j_groups[[i]]]
            values <- paste0(values[-which(values == "N")], collapse = "")
            temp <- temp[startsWith(temp$codon, values), ]
          }
        } else{
          values <- paste0(ref_j[j_groups[[i]]], collapse = "")
          temp <- temp[endsWith(temp$codon, values), ]
        }
      } else{
        if(length(j_groups[[i]]) == 3){
          if(sum("N" %in% ref_j[j_groups[[i]]]) == 0){
            if(alakazam::translateDNA(paste0(ref_j[j_groups[[i]]], collapse = "")) != "*"){
              temp <- temp[temp$codon == paste0(ref_j[j_groups[[i]]], collapse = "") |
                             (alakazam::translateDNA(temp$codon) %in% c("F", "W")),]
            } else{
              temp <- temp[alakazam::translateDNA(temp$codon) %in% c("F", "W"),]
            }
          } else{
            values <- ref_j[j_groups[[i]]]
            values <- paste0(values[-which(values == "N")], collapse = "")
            temp <- temp[(startsWith(temp$codon, values)) |
                           (alakazam::translateDNA(temp$codon) %in% c("F", "W")), ]
          }
        } else{
          values <- paste0(ref_j[j_groups[[i]]], collapse = "")
          temp <- temp[(endsWith(temp$codon, values)) |
                         (alakazam::translateDNA(temp$codon) %in% c("F", "W")),]
        }
      }
    } else{
      if(length(j_groups[[i]]) == 3){
        if(sum("N" %in% ref_j[j_groups[[i]]]) == 0){
          temp <- temp[temp$codon == paste0(ref_j[j_groups[[i]]], collapse = ""),]
        } else{
          values <- ref_j[j_groups[[i]]]
          values <- paste0(values[-which(values == "N")], collapse = "")
          temp <- temp[startsWith(temp$codon, values), ]
        }
      } else{
        values <- paste0(ref_j[j_groups[[i]]], collapse = "")
        temp <- temp[endsWith(temp$codon, values), ]
      }
    }
    return(temp)
  }))
  j_df <- j_df[, !names(j_df) %in% "new_site"]
  junc_df <- tree_df[!tree_df$site %in% c(unique(v_df$site), unique(j_df$site)),]
  write.table(tree_df, file.path(subDir, paste0("original_", clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")), 
              quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  tree_df <- rbind(v_df, junc_df, j_df)
  tree_seq <- paste0(unlist(lapply(unique(tree_df$site), function(x){
    sub <- tree_df[tree_df$site == x,]
    sub_indx <- which(sub$value == max(sub$value))[1]
    sub$codon[sub_indx]
  })), collapse = "")
  tree_df <- tree_df[, !names(tree_df) == "value"]
  if(sub$data[[1]]@phylo_seq != "hlsequence"){
    write.table(tree_df, file.path(subDir, paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")), 
                quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  } else{
    if(chain == "H"){
      write.table(tree_df, file.path(subDir, "heavy_table.txt"), 
                  quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
    } else{
      write.table(tree_df, file.path(subDir, "light_table.txt"), 
                  quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
    }
  }
  # update the starting values to reflect these changes 
  v <- substring(tree_seq, 1, nchar(v))
  test_junc <- substring(tree_seq, nchar(v) + 1, nchar(v) + nchar(mrcacdr3))
  cdr3_c <- substring(mrcacdr3, 1, 3)
  cdr3_fw <- substring(mrcacdr3, nchar(mrcacdr3)-2, nchar(mrcacdr3))
  mrcacdr3 <- paste0(cdr3_c, substring(test_junc, 4, nchar(test_junc)-3),
                     cdr3_fw)
  if(pad_length == 0){
    j <- substring(tree_seq, nchar(v) + nchar(mrcacdr3) + 1, nchar(tree_seq))
  } else{
    j <- substring(tree_seq, nchar(v) + nchar(mrcacdr3) + 1, nchar(tree_seq)- pad_length)
    j <- paste0(j, paste(rep("N", pad_length), collapse = ""))
  }
  temp <- data.frame(v = v, j = j, cdr3 = mrcacdr3)
  return(temp)
}

# Prepares clones for UCA inference. 
#
# \code{processCloneGermline} Exports two text files that are used as inputs for 
#                             the UCA script.
# @param clone_ids  A clone id for the clones object. This is only used in parallel
# @param clones     A clones object from \link{formatClones}
# @param data       The airr-table associated with the clones object
# @param dir        The directory where data should be saved to 
# @param id         The run id
# @param quiet      How much noise to print out
# @param chain      HL or H?
# @param check_genes Check if the inferred V/J lengths go into the inferred cdr3 region and adjust accordingly. 
# @param references IMGT references read in using \link{readIMGT}
#
processCloneGermline <- function(clone_ids, clones, data, dir, build, id,
                                 quiet = 0, chain = "H", check_genes = FALSE,
                                 references = NULL){
  sub <- dplyr::filter(clones, !!rlang::sym("clone_id") == clone_ids)
  subDir <- file.path(dir, paste0(id, "_",clone_ids))
  if(!dir.exists(subDir)){
    dir.create(subDir, recursive = T)
  }
  saveRDS(sub, file.path(subDir, "clone.rds"))
  
  if(sub$data[[1]]@phylo_seq == "hlsequence"){
    test_hl <- paste0(strsplit(sub$data[[1]]@hlgermline, "")[[1]][(nchar(sub$data[[1]]@germline) + 1):
                                                             nchar(sub$data[[1]]@hlgermline)], collapse = "")
    if(nchar(test_hl) > nchar(sub$data[[1]]@lgermline)){
      sub$data[[1]]@lgermline <- test_hl
    }
  }
  
  imgt_germline <- sub$data[[1]]@germline
  if(sub$data[[1]]@phylo_seq == "lsequence"){
    imgt_germline <- sub$data[[1]]@lgermline
  }
  r <- sub$data[[1]]@region
  if(sub$data[[1]]@phylo_seq == "sequence"){
    cdr3_index <- (min(which(r == "cdr3")) - 3):(max(which(r == "cdr3")) + 3)
  } else if(sub$data[[1]]@phylo_seq == "hlsequence"){
    heavy_r <- r[1:nchar(sub$data[[1]]@germline)]
    light_r <- r[(nchar(sub$data[[1]]@germline) + 1): length(r)]
    cdr3_index <- (min(which(heavy_r == "cdr3")) - 3):(max(which(heavy_r == "cdr3")) + 3)
  } else if(sub$data[[1]]@phylo_seq == "lsequence"){
    cdr3_index <- (min(which(r == "cdr3")) - 3):(max(which(r == "cdr3")) + 3)
  }
  if(build == "igphyml"){
    tree_df <- suppressWarnings(read.table(file = file.path(dir, "sample", "sample_recon_sample",
                                           paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")), 
                                           header = F, sep = "\t"))
    file.copy(file.path(dir, "sample", "sample_recon_sample",
                        paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")), 
              file.path(subDir, paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")))
  } else if(build == "pml"){
    tree_df <- suppressWarnings(read.table(file = file.path(subDir, "codon_table.txt"), 
                                           header = F, sep = "\t"))
    file.rename(file.path(subDir, "codon_table.txt"), 
                file.path(subDir, paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")))
  }

  colnames(tree_df) = c("site", "codon", "partial_likelihood", "log_likelihood_site",
                        "upper_partial_log_likelihood", "upper_partial_likelihood", "equilibrium")
  tree_df$value <- tree_df$partial_likelihood + log(tree_df$equilibrium)
  
  tree_seq <- paste0(unlist(lapply(unique(tree_df$site), function(x){
    sub <- tree_df[tree_df$site == x,]
    sub_indx <- which(sub$value == max(sub$value))[1]
    sub$codon[sub_indx]
  })), collapse = "")
  
  if(sub$data[[1]]@phylo_seq == "hlsequence"){
    nsite_heavy <- nchar(sub$data[[1]]@germline)/3
    tree_df_light <- tree_df[tree_df$site >= nsite_heavy,]
    tree_df_light$site <- tree_df_light$site - min(tree_df_light$site)
    tree_df <- tree_df[tree_df$site < nsite_heavy,]
    mrca <- substring(tree_seq, 1, nchar(sub$data[[1]]@germline))
    mrcacdr3 <- paste0(strsplit(mrca, "")[[1]][cdr3_index], collapse = "")
    mrca_light <- substring(tree_seq, nchar(sub$data[[1]]@germline) + 1, nchar(tree_seq))
  }else if(sub$data[[1]]@phylo_seq == "sequence"){
    mrcacdr3 <- paste0(strsplit(tree_seq, "")[[1]][cdr3_index], collapse = "")
  } else{
    mrcacdr3 <- paste0(strsplit(tree_seq, "")[[1]][cdr3_index], collapse = "")
  }
  
  test_cdr3 <- strsplit(alakazam::translateDNA(mrcacdr3), "")[[1]]
  if(test_cdr3[1] != "C" || !test_cdr3[length(test_cdr3)] %in% c("F", "W")){
    if(sub$data[[1]]@phylo_seq == "sequence"){
      groupedList <- split(1:length(r), ceiling(seq_along(r) / 3))
    } else if(sub$data[[1]]@phylo_seq == "hlsequence"){
      groupedList <- split(1:length(heavy_r), ceiling(seq_along(heavy_r) / 3))
    } else if(sub$data[[1]]@phylo_seq == "lsequence"){
      groupedList <- split(1:length(r), ceiling(seq_along(r) / 3))
    }
    if(test_cdr3[1] != "C"){
      codon_site <- which(sapply(groupedList, function(group) min(cdr3_index) %in% group))
      sub_tree_df <- dplyr::filter(tree_df, !!rlang::sym("site") == codon_site - 1)
      sub_tree_df$aa <- alakazam::translateDNA(sub_tree_df$codon)
      sub_tree_df <- dplyr::filter(sub_tree_df, !!rlang::sym("aa") == "C")
      sub_tree_df$value <- sub_tree_df$partial_likelihood + log(sub_tree_df$equilibrium)
      value <- sub_tree_df$codon[sub_tree_df$value == max(sub_tree_df$value)]
      mrcacdr3 <- paste0(value[1], substring(mrcacdr3, 4, nchar(mrcacdr3)))
    }
    if(!test_cdr3[length(test_cdr3)] %in% c("F", "W")){
      codon_site <- which(sapply(groupedList, function(group) max(cdr3_index) %in% group))
      sub_tree_df <- dplyr::filter(tree_df, !!rlang::sym("site") == codon_site - 1)
      sub_tree_df$aa <- alakazam::translateDNA(sub_tree_df$codon)
      sub_tree_df <- dplyr::filter(sub_tree_df, !!rlang::sym("aa") %in% c("W", "F"))
      sub_tree_df$value <- sub_tree_df$partial_likelihood + log(sub_tree_df$equilibrium)
      value <- sub_tree_df$codon[sub_tree_df$value == max(sub_tree_df$value)]
      mrcacdr3 <- paste0(substring(mrcacdr3, 1, nchar(mrcacdr3)-3), value[1])
    }
  }
  
  if(sub$data[[1]]@phylo_seq == "sequence"){
    v_len <- min(cdr3_index)-1
    v <- substring(imgt_germline, 1, v_len)
    j_start <- nchar(paste0(v, mrcacdr3, collapse = "")) +1
    j <- substring(imgt_germline, j_start, nchar(imgt_germline))
    j_len <- nchar(j)
  } else if(sub$data[[1]]@phylo_seq == "hlsequence"){
    v_len <- min(cdr3_index)-1
    v <- substring(imgt_germline, 1, v_len)
    j_start <- nchar(paste0(v, mrcacdr3, collapse = "")) +1
    j <- substring(sub$data[[1]]@germline, j_start, nchar(sub$data[[1]]@germline))
    j_len <- nchar(j)
  } else if(sub$data[[1]]@phylo_seq == "lsequence"){
    v_len <- min(cdr3_index)-1
    v <- substring(imgt_germline, 1, v_len)
    j_start <- nchar(paste0(v, mrcacdr3, collapse = "")) +1
    j <- substring(imgt_germline, j_start, nchar(imgt_germline))
    j_len <- nchar(j)
  }
  if(quiet > 0){
    print(paste("sucessfully obtained most likely junction for", clone_ids))
  }

  if(sub$data[[1]]@phylo_seq == "hlsequence"){
    v_light <- substring(sub$data[[1]]@lgermline, 1, sum(light_r %in% c("cdr1", "cdr2", "fwr1", "fwr2", "fwr3")))
    light_cdr3 <- substring(sub$data[[1]]@lgermline, nchar(v_light) + 1, nchar(v_light) + sum(light_r == "cdr3"))
    j_light <- substring(sub$data[[1]]@lgermline, nchar(v_light) + nchar(light_cdr3) + 1, nchar(sub$data[[1]]@lgermline))
    
    last_v_codon <- substring(v_light, nchar(v_light)-2, nchar(v_light))
    first_j_codon <- substring(j_light, 1, 3)
    light_cdr3 <- paste0(last_v_codon, light_cdr3, first_j_codon)
    v_light <- substring(v_light, 1, nchar(v_light)-3)
    j_light <- substring(j_light, 4, nchar(j_light))
    
    light_cdr3_test <- strsplit(alakazam::translateDNA(light_cdr3), "")[[1]]
    
    if(light_cdr3_test[1] != "C" || !light_cdr3_test[length(light_cdr3_test)] %in% c("F", "W")){
      if(light_cdr3_test[1] != "C"){
        codon_site <- nchar(v_light)/3 + 1
        sub_tree_df <- dplyr::filter(tree_df_light, !!rlang::sym("site") == codon_site)
        sub_tree_df$aa <- alakazam::translateDNA(sub_tree_df$codon)
        sub_tree_df <- dplyr::filter(sub_tree_df, !!rlang::sym("aa") == "C")
        sub_tree_df$value <- sub_tree_df$partial_likelihood + log(sub_tree_df$equilibrium)
        value <- sub_tree_df$codon[sub_tree_df$value == max(sub_tree_df$value)]
        light_cdr3 <- paste0(value[1], substring(light_cdr3, 4, nchar(light_cdr3)))
      }
      if(!light_cdr3_test[length(light_cdr3_test)] %in% c("F", "W")){
        codon_site <- nchar(paste0(v_light, light_cdr3))/3 + 1
        sub_tree_df <- dplyr::filter(tree_df_light, !!rlang::sym("site") == codon_site)
        sub_tree_df$aa <- alakazam::translateDNA(sub_tree_df$codon)
        sub_tree_df <- dplyr::filter(sub_tree_df, !!rlang::sym("aa") %in% c("W", "F"))
        sub_tree_df$value <- sub_tree_df$partial_likelihood + log(sub_tree_df$equilibrium)
        value <- sub_tree_df$codon[sub_tree_df$value == max(sub_tree_df$value)]
        light_cdr3 <- paste0(substring(light_cdr3, 1, nchar(light_cdr3)-3), value[1])
      }
    }
  }
  saveRDS(sub, file.path(subDir, "clone.rds"))
  if(check_genes){
    if(sub$data[[1]]@phylo_seq == "hlsequence"){
      heavy_vals <- checkGenesUCA(sub = sub, data = data, v = v, mrcacdr3 = mrcacdr3,
                                  j = j, references = references, tree_df = tree_df,
                                  subDir = subDir, clone_ids = clone_ids, chain = "H",
                                  r = heavy_r)
      light_vals <- checkGenesUCA(sub = sub, data = data, v = v_light, mrcacdr3 = light_cdr3,
                                  j = j_light, references = references, tree_df = tree_df_light,
                                  subDir = subDir, clone_ids = clone_ids, chain = "L",
                                  r = light_r)
      v <- heavy_vals$v
      mrcacdr3 <- heavy_vals$cdr3
      j <- heavy_vals$j
      v_light <- light_vals$v
      light_cdr3 <- light_vals$cdr3
      j_light <- light_vals$j
    } else if(sub$data[[1]]@phylo_seq == "sequence"){
      heavy_vals <- checkGenesUCA(sub = sub, data = data, v = v, mrcacdr3 = mrcacdr3,
                                  j = j, references = references, tree_df = tree_df,
                                  subDir = subDir, clone_ids = clone_ids, chain = "H",
                                  r = r)
      v <- heavy_vals$v
      mrcacdr3 <- heavy_vals$cdr3
      j <- heavy_vals$j
    } else if(sub$data[[1]]@phylo_seq == "lsequence"){
      light_vals <- checkGenesUCA(sub = sub, data = data, v = v, mrcacdr3 = mrcacdr3,
                                  j = j, references = references, tree_df = tree_df,
                                  subDir = subDir, clone_ids = clone_ids, chain = "L",
                                  r = r)
      v <- light_vals$v
      mrcacdr3 <- light_vals$cdr3
      j <- light_vals$j
    }
  }
  # put it all together 
  v_cdr3 <- paste0(v, paste0(mrcacdr3, collapse = ""), collapse = "")
  starting_germ <- paste0(v_cdr3, j, collapse = "")
  file_path_germline <- file.path(subDir, paste("olga_testing_germline.txt"))
  file_path_junction_position <- file.path(subDir, paste("olga_junction_positions.txt"))
  writeLines(paste0(starting_germ, collapse = ""), con = file_path_germline)
  writeLines(paste(min(cdr3_index)-1, max(cdr3_index)), con = file_path_junction_position)
  if(sub$data[[1]]@phylo_seq == "hlsequence"){
    v_cdr3 <- paste0(v_light, light_cdr3, collapse = "")
    starting_germ <- paste0(v_cdr3, j_light, collapse = "")
    file_path_germline <- file.path(subDir, paste("olga_testing_germline_light.txt"))
    file_path_junction_position <- file.path(subDir, paste("olga_junction_positions_light.txt"))
    writeLines(paste0(starting_germ, collapse = ""), con = file_path_germline)
    writeLines(paste(nchar(v_light), nchar(v_cdr3)), con = file_path_junction_position)
    if(!check_genes){
      tree_df <- tree_df[, !names(tree_df) %in% c("value")]
      write.table(tree_df, file.path(subDir, "heavy_table.txt"), quote = FALSE,
                  sep = "\t", col.names = FALSE, row.names = FALSE)
      tree_df_light <- tree_df_light[, !names(tree_df_light) %in% c("value")]
      write.table(tree_df_light, file.path(subDir, "light_table.txt"), quote = FALSE,
                  sep = "\t", col.names = FALSE, row.names = FALSE)
    }

  }
  return(sub)
}

# Runs clones through a UCA inference. 
#
# \code{callOlga} Performs UCA inference and exports the UCA, some data about the UCA,
#                 and UCA likelihoods
# @param clones     The airrClones object
# @param dir        The directory where data should be saved to 
# @param uca_script The file path to the UCA python script
# @param python     The call used to launch python from command line
# @param max_iters  The maximum number of iterations to run before ending
# @param nproc      The number of cores to use 
# @param id         The run id
# @param model_folder The file path to the model parameters for IGH provide by OLGA
# @param model_folder_igk The file path to the model parameters for IGK provide by OLGA
# @param model_folder_igl The file path to the model parameters for IGL provide by OLGA
# @param quiet      Amount of noise to print out
# @param search     Search the codon or nt space
# @param clone      The column name for the clone identifier 
#

callOlga <- function(clones, dir, uca_script, python, max_iters, nproc, id, model_folder,
                     model_folder_igk, model_folder_igl, quiet, search, clone){
  clone_ids <- paste0(unlist(lapply(clones[[clone]], function(z){
    value <- z
    if(clones$data[[which(clones[[clone]] == z)]]@phylo_seq == "hlsequence"){
      value <- append(value, z)
    }
    return(value)
  })), collapse = ",")
  starting_germlines <- paste0(unlist(lapply(clones[[clone]], function(z){
    value <- path.expand(file.path(dir, paste0(id, "_", z), "olga_testing_germline.txt"))
    if(clones$data[[which(clones[[clone]] == z)]]@phylo_seq == "hlsequence"){
      value <- append(value, path.expand(file.path(dir, paste0(id, "_", z), "olga_testing_germline_light.txt")))
    }
    return(value)
  })), collapse = ",")
  junction_location <- paste0(unlist(lapply(clones[[clone]], function(z){
    value <- path.expand(file.path(dir, paste0(id, "_", z), "olga_junction_positions.txt"))
    if(clones$data[[which(clones[[clone]] == z)]]@phylo_seq == "hlsequence"){
      value <- append(value, path.expand(file.path(dir, paste0(id, "_", z), "olga_junction_positions_light.txt")))
    }
    return(value)
  })), collapse = ",")
  tree_tables <- paste0(unlist(lapply(clones[[clone]], function(z){
    if(clones$data[[which(clones[[clone]] == z)]]@phylo_seq == "sequence"){
      value <- path.expand(file.path(dir, paste0(id, "_", z), 
                                     paste0(z, ".fasta_igphyml_rootprobs_hlp.txt")))
    } else if(clones$data[[which(clones[[clone]] == z)]]@phylo_seq == "hlsequence"){
      value <- path.expand(file.path(dir, paste0(id, "_", z), 
                                     "heavy_table.txt"))
      value <- append(value, path.expand(file.path(dir, paste0(id, "_", z),
                                                   "light_table.txt")))
    } else if(clones$data[[which(clones[[clone]] == z)]]@phylo_seq == "lsequence"){
      value <- path.expand(file.path(dir, paste0(id, "_", z), 
                                     paste0(z, ".fasta_igphyml_rootprobs_hlp.txt")))
    }
    return(value)
  })), collapse = ",")
  chains <- paste0(unlist(lapply(clones[[clone]], function(z){
    loci <- strsplit(clones$locus[which(clones[[clone]] == z)], ",")[[1]]
  })), collapse = ",")
  
  args <- c(
    "--clone_ids", clone_ids, 
    "--directory", path.expand(dir), 
    "--max_iters", max_iters,
    "--nproc", nproc, 
    "--id", id, 
    "--model_folder", path.expand(model_folder),
    "--model_folder_igk", ifelse(is.null(model_folder_igk), "NULL", path.expand(model_folder_igk)),
    "--model_folder_igl", ifelse(is.null(model_folder_igl), "NULL", path.expand(model_folder_igl)),
    "--quiet", quiet,
    "--starting_germlines", starting_germlines,
    "--junction_locations", junction_location,
    "--tree_tables", tree_tables,
    "--chains", chains,
    "--search", search
  )
  cmd <- paste(
    shQuote(python), 
    shQuote(path.expand(uca_script)), 
    paste(shQuote(args), collapse=" ")
  )
  olga_check <- tryCatch(system(cmd), error=function(e)e)
  if("error" %in% class(olga_check)){
      stop("there was an error running the get_UCA script. This is likely due to 
           not having the required python packages installed for the python version found
           at", path.expand(Sys.which(python)))
  }
}

# Adds the UCA to the clones object 
#
# \code{updateClone} Adds the UCA to the data frame within the clones object
# @param clones     The clones object
# @param dir        The directory where data should be saved to 
# @param id         The run id
# @param nproc      The number of cores to use 
#

updateClone <- function(clones, dir, id, nproc = 1){
  updated_clones <- do.call(rbind, parallel::mclapply(1:nrow(clones), function(x){
    clone <- clones[x,]
    if(clone$data[[1]]@phylo_seq == "lsequence"){
      uca <- read.table(file.path(dir, paste0(id, "_", clone$clone_id),
                                  "UCA_light.txt"), sep = "\t")[[1]]
    } else{
      uca <- read.table(file.path(dir, paste0(id, "_", clone$clone_id), "UCA.txt"), sep = "\t")[[1]]
      if(clone$data[[1]]@phylo_seq == "hlsequence"){
        uca <- paste0(uca, read.table(file.path(dir, paste0(id, "_", clone$clone_id),
                                                "UCA_light.txt"), sep = "\t")[[1]])
      }
    }
    clone$UCA <- uca
    clone$AA_UCA <- alakazam::translateDNA(uca)
    germline_node <- ape::getMRCA(clone$trees[[1]], clone$trees[[1]]$tip.label)
    clone$trees[[1]]$nodes[[germline_node]]$sequence <- uca
    if(clone$data[[1]]@phylo_seq == "hlsequence"){
      UCA_heavy <- getNodeSeq(clone, node = germline_node, tree = clone$trees[[1]])[1]
      UCA_light <- getNodeSeq(clone, node = germline_node, tree = clone$trees[[1]])[2]
      clone$UCA_IMGT <- I(list(c(UCA_heavy, UCA_light)))
    } else if(clone$data[[1]]@phylo_seq == "sequence"){
      clone$UCA_IMGT <- getNodeSeq(clone, node = germline_node, tree = clone$trees[[1]])[1]
    } else if(clone$data[[1]]@phylo_seq == "lsequence"){
      clone$UCA_IMGT <- getNodeSeq(clone, node = germline_node, tree = clone$trees[[1]])[1]
    }
    return(clone)
  }, mc.cores = nproc))
  return(updated_clones)
}

#' \link{getTreesAndUCAs} Construct trees and infer the UCA
#' @param clones        AIRR-table containing sequences \link{formatClones}
#' @param data          The AIRR-table that was used to make the clones object.
#' @param dir           The file path of the directory of where data is saved. NULL is default.
#' @param build         Name of the tree building method
#' @param exec          File path to the tree building executable
#' @param model_folder  The file path to the OLGA default model files for heavy chains
#' @param model_folder_igk  The file path to the OLGA default model files for IGK
#' @param model_folder_igl  The file path to the OLGA default model files for IGL
#' @param uca_script    The file path to the UCA python script
#' @param python        Specify the python call for your system. This is the call
#'                      on command line that issues the python you want to use. 
#'                      "python" by default. 
#' @param id            The run ID
#' @param max_iters     The maximum number of iterations to run before ending
#' @param nproc         The number of cores to use 
#' @param rm_temp       Remove the generated files?
#' @param chain         Set to HL to use both heavy and light chain sequences
#' @param quiet         Amount of noise to print out
#' @param references    Reference genes. See \link{readIMGT}
#' @param clone         The name of the clone id column used in \link{formatClones}
#' @param cell          The name of the cell id in the AIRR table used to generate \link{formatClones}
#' @param heavy         The name of the heavy chain locus. Default is IGH. 
#' @param sampling_method How to subsample. Methods include 'random', 'lm' (least mutated),
#'                      and 'ratio' (a weighted sampling with heavier weights
#'                      towards the least mutated sequences). The later two methods
#'                      require 'mu_freq' to be passed as a trait when running
#'                      \link{formatClones}
#' @param subsample_size The amount that the clone should be sampled down to. Default is NA. Use NA if you do not wish to subsample -- in testing
#' @param search        Search codon or nt space
#' @param check_genes   Check if the inferred V/J lengths go into the inferred cdr3 region and adjust accordingly. 
#' @param igblast       File path to igblast 
#' @param igblast_database The file path to the database setup for igblast 
#' @param ref_path      The file path to your references parent folder.
#' @param organism      The type of organism to align to if using igblast. 
#' @param ...           Additional arguments passed to various other functions like \link{getTrees} and \link{buildGermline}
#' @return An \code{airrClone} object with trees and the inferred UCA
#' @details Return object adds/edits following columns:
#' \itemize{
#'   \item  \code{trees}:  The phylogenies associated with each clone
#'   \item  \code{UCA}:    The inferred UCA
#' }
#' @seealso \link{getTrees} 
#' @export
getTreesAndUCAs <- function(clones, data, dir = NULL, build, exec,  model_folder,
                           model_folder_igk = NULL, model_folder_igl = NULL, 
                           uca_script, python = "python", id = "sample", 
                           max_iters = 100, nproc = 1, rm_temp = TRUE, quiet = 0,
                           chain = "H", references = NULL, clone = "clone_id",
                           cell = "cell_id", heavy = "IGH", sampling_method = 'random',
                           subsample_size = NA, search = "codon", check_genes = TRUE,
                           igblast = NULL, igblast_database = NULL, ref_path = NULL, 
                           organism = 'human', ...){
  if(!is.null(dir)){
    dir <- path.expand(dir)
    dir <- file.path(dir, paste0("all_", id))
    if(!dir.exists(dir)){
      dir.create(dir, recursive = TRUE)
    }
  }else{
    dir <- alakazam::makeTempDir(id)
  }
  
  if(rm_temp){
    rm_dir = dir
  } else{
    rm_dir = NULL
  }
  
  if(file.access(path.expand(Sys.which(python)), mode = 1) == -1){
    stop("The python executable found at", path.expand(Sys.which(python)), " cannot be executed.")
  }
  
  if(!file.exists(path.expand(uca_script))){
    stop("The file", path.expand(uca_script), " cannot be executed.")
  }
  
  if(check_genes & is.null(references)){
    stop("check_genes cannot be run without references. Pass a references object using references call.",
         "References need to be read in using dowser::readIMGT()")
  }
  
  if(!is.null(igblast) & any(sapply(list(igblast_database, ref_path), is.null))){
    stop('if you are planning on updating the germline references using igblast, the following varaibles cannot be NULL:',
         "\n",
         'igblast: the file path to igblast', "\n",
         'igblast_database: the file path to the databases set up for igblast',"\n",
         'ref_path: the path to the imgt parent folder')
  }
  if(!is.na(subsample_size)){
    if(!is.numeric(subsample_size)){
      stop("subsample_size must be a numeric")
    }
    if(sampling_method == "random"){
      clones <- sampleClones(clones, size = subsample_size)
    } else{
      if(!"mu_freq" %in% colnames(clones$data[[1]]@data)){
        stop('Mutation frequency calculations are required for this subsampling method.',
             ' Please run your data through shazam::observedMutations() and then rerun',
             ' formatClones and getTreesAndUCAs.')
      }
      if(sampling_method == "lm"){
        clones <- do.call(rbind, parallel::mclapply(1:nrow(clones), function(x){
          sub <- clones[x,]
          if(!subsample_size > nrow(sub$data[[1]]@data)){
            temp_subsample_size <- subsample_size
          } else{
            temp_subsample_size <- nrow(sub$data[[1]]@data)
          }
          sub$data[[1]]@data$weight <- NA
          temp_df <- sub$data[[1]]@data[order(sub$data[[1]]@data$mu_freq, decreasing = FALSE),]
          sample_names <- temp_df$sequence_id[1:temp_subsample_size]
          sub$data[[1]]@data$weight[sub$data[[1]]@data$sequence_id %in% sample_names] <- 100
          sub$data[[1]]@data$weight[!sub$data[[1]]@data$sequence_id %in% sample_names] <- 0
          return(sub)
        }, mc.cores = nproc))
        clones <- sampleClones(clones, subsample_size, weight = "weight")
      } else if(sampling_method == "ratio"){
        clones <- do.call(rbind, parallel::mclapply(1:nrow(clones), function(x){
          sub <- clones[x,]
          sub$data[[1]]@data$weight <- 1/(sub$data[[1]]@data$mu_freq + 1e-10)
          return(sub)
        }, mc.cores = nproc))
        clones <- sampleClones(clones, subsample_size, weight = "weight")
      } else{
        stop('sampling_method:', sampling_method, ' not recognized')
      }
    }
    # if subsampling to 1 sequence
    if(subsample_size == 1){
      cells <- unlist(lapply(clones$data, function(x) x@data$sequence_id))
      if (!is.na(cell) & cell %in% colnames(data)) {
        filtered <- data[data$sequence_id %in% cells, ]
        cells <- data$sequence_id[data[[cell]] %in% filtered[[cell]]]
      }
      sub_data <- data[data$sequence_id %in% cells,]
      clones <- formatClones(sub_data, nproc = nproc, 
                             filterstop = TRUE, chain = chain, minseq = 1, 
                             dup_singles = T, traits = "mu_freq")
    }
    saveRDS(clones, file = file.path(dir, "clones.rds"))
  }
  if(quiet > 0){
    print("constructing trees")
  }
  if(build == "igphyml"){
    clones <- getTrees(clones, build = build, exec = exec, rm_temp = FALSE, dir = dir,
                    asrp = TRUE, chain = chain, nproc = nproc, ...)
  } else if(build == "pml"){
    clones <- getTrees(clones, build = build, rm_temp = FALSE, dir = dir, asrp = TRUE,
                       nproc = nproc, ...)
  } else{
    stop("the tree bulding method ", build, "is not supported")
  }
  saveRDS(clones, file.path(dir, "clones.rds"))
  if(!is.null(igblast)){
    if(quiet > 0){
      print("updating cuts based on MRCA")
    }
    data <- updateAIRRGerm(airr_data = data, clones = clones, igblast = igblast, 
                           igblast_database = igblast_database, references = ref_path, 
                           organism = organism, locus = 'Ig', outdir = dir, 
                           nproc = nproc)
    airr::write_rearrangement(data, file = file.path(dir, "updated_data.tsv"))
  }
  if(quiet > 0){
    print("preparing the clones for UCA analysis")
  }
  clones <- do.call(rbind, invisible(parallel::mclapply(clones[[clone]], function(x){
    processCloneGermline(clone_ids = x, clones = clones, data = data, dir = dir,
                         build = build, id = id, quiet = quiet, 
                         chain = chain, check_genes = check_genes, 
                         references = references)
  }, mc.cores = nproc)))
  # run the UCA
  if(quiet > 0){
    print("running UCA analysis")
  }
  callOlga(clones = clones, dir = dir, model_folder = model_folder,
           model_folder_igk = model_folder_igk, model_folder_igl = model_folder_igl,
           uca_script = uca_script, python = python, max_iters = max_iters, nproc = nproc,
           id = id, quiet = quiet, searc = search, clone = clone)

  if(quiet > 0){
    print("updating clones")
  }
  clones <- updateClone(clones, dir, id, nproc)

  unlink(rm_dir,recursive=TRUE)
  return(clones)
}
