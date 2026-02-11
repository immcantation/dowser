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
assignGenes <- function(file, igblast, refs, igdata = NULL, organism = "human",
    domain_system = "imgt", outfile = NULL, nproc = 1, db_prefix="imgt",
    locus = "Ig", set_igdata=TRUE, return=TRUE, verbose=TRUE){
  
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
# param clone the name of the proper clone_id column to use
# param v     the v string from resolving the germline previously
# param j     the j string from resolving the germline previously 
# param v_light the v string from resolving the germline previously. Only for HL trees
# param j_light the j string from resolving the germline previously. Only for HL trees 
# return directory of downloaded IMGT BCR and TCR reference files
# 
updateAIRRGerm <- function(airr_data, clones, igblast, igblast_database, 
                           references, igdata = NULL, organism = "human", 
                           locus = "Ig", outdir, nproc = 1, clone = "clone_id", 
                           v = NULL, j = NULL, v_light = NULL, j_light = NULL,...){
  if(is.null(igdata)){
    igdata <- igblast
  }
  mrcas <- do.call(rbind, parallel::mclapply(1:nrow(clones), function(x){
    if(clones$data[[x]]@phylo_seq != "hlsequence"){
      value <- getMRCASeq(clones[x,])
      if(!is.null(v) & !is.null(j)){
        mrca_junc <- substring(value, nchar(v)+1, nchar(value) - nchar(j))
        value <- paste0(v, mrca_junc, j)
      }
      if(clones$data[[x]]@phylo_seq == "sequence"){
        seq_id <- paste0(clones$clone_id[x], "_heavy")
        og_germ <- airr_data$germline_alignment[airr_data[[clone]] == clones$clone_id[x] &
                                                  airr_data$locus == "IGH"][1]
      } else{
        seq_id <- paste0(clones$clone_id[x], "_light")
        og_germ <- airr_data$germline_alignment[airr_data[[clone]] == clones$clone_id[x] &
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
      if(!is.null(v) & !is.null(j)){
        mrca_junc <- substring(hvalue, nchar(v)+1, nchar(hvalue) - nchar(j))
        hvalue <- paste0(v, mrca_junc, j)
      }
      lvalue <- substring(value, length(hnumbers)+1, nchar(value))
      if(!is.null(v_light) & !is.null(j_light)){
        mrca_junc <- substring(lvalue, nchar(v_light)+1, nchar(lvalue) - nchar(j_light))
        lvalue <- paste0(v_light, mrca_junc, j_light)
      }
      seq_id <- c(paste0(clones$clone_id[x], "_heavy"), paste0(clones$clone_id[x], "_light"))

      og_germ <- airr_data$germline_alignment[airr_data[[clone]] == clones$clone_id[x] &
                                                airr_data$locus == "IGH"][1]
      gaps <- dplyr::setdiff(1:max(hnumbers), hnumbers)
      og_germ <- paste0(strsplit(og_germ, "")[[1]][-gaps], collapse = "")
      og_lgerm <- airr_data$germline_alignment[airr_data[[clone]] == clones$clone_id[x] &
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
                         igdata = igdata, 
                         outfile = file.path(outdir, "igblast.tsv"),
                         nproc = nproc, ...)
  ig_data <- addGaps(ig_data, gapdb=references, organism=organism, locus=locus, ...)
  airr::write_rearrangement(ig_data, file.path(outdir, "igblast.tsv"))
  ig_data[[clone]] <- substring(ig_data$sequence_id, 1, nchar(ig_data$sequence_id) - 6)

  # get the length values 
  ig_data$v_germline_length <- as.numeric(ig_data$v_germline_end - ig_data$v_germline_start + 1)
  ig_data$d_germline_length <- as.numeric(ig_data$d_germline_end - ig_data$d_germline_start + 1)
  ig_data$j_germline_length <- as.numeric(ig_data$j_germline_end - ig_data$j_germline_start + 1)
  
  if(sum(is.na(ig_data$j_germline_length)) != 0){
    bad_clones <- ig_data[[clone]][which(is.na(ig_data$j_germline_length))]
    tmp <- airr_data[airr_data[[clone]] %in% bad_clones,]
    ig_data <- ig_data[!ig_data[[clone]] %in% bad_clones,]
    if(nrow(ig_data) == 0){
      # return the the orignal data as the MRCA failed 
      return(airr_data)
    }
  }
  
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
  airr::write_rearrangement(ig_data, file.path(outdir, "igblast.tsv"))
  cols_to_update <- c("germline_alignment", "germline_alignment_d_mask", "v_germline_start", "v_germline_end",
                      "d_germline_start", "d_germline_end", "j_germline_start", "j_germline_end",
                      "np1_length", "np2_length", "v_germline_length", "d_germline_length", "j_germline_length",
                      "v_call", "d_call", "j_call")
  airr_data <- merge(airr_data, ig_data[, c(clone, "locus", cols_to_update)], 
                     by = c(clone, "locus"), all.x = TRUE, suffixes = c("", ".ig"))
  for (col in cols_to_update) {
    ig_col <- paste0(col, ".ig")
    if(ig_col %in% names(airr_data)){
      tmp <- airr_data[[col]]
      airr_data[[col]] <- airr_data[[ig_col]]
      airr_data[[ig_col]] <- tmp
    }
  }
  airr_data <- airr_data[!is.na(airr_data$v_call),]
  airr::write_rearrangement(airr_data, file.path(outdir, "updated_airr_data.tsv"))
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
# clone is the name of the clone id column to use
# seq_id is the name of the sequence id column to use
checkGenesUCA <- function(sub, data, v, mrcacdr3, j, references, tree_df, subDir,
                          clone_ids, r, chain = "H", clone = "clone_id", 
                          seq_id = "sequence_id"){
  if(chain == "H"){
    cons <- findConsensus(receptors = data[data[[clone]] == sub$clone_id & data$locus == "IGH",],
                          clone_id = sub$clone_id)
  } else{
    cons <- findConsensus(receptors = data[data[[clone]] == sub$clone_id & data$locus != "IGH",],
                          clone_id = sub$clone_id)
  }
  cons <- data[data[[seq_id]] == cons$cons_id,]
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
    cons$j_germline_length <- as.numeric(cons$j_germline_length) + ig_diff
    cons$j_germline_end <- as.numeric(cons$j_germline_end) + ig_diff
  }
  ref_j <- substring(ref_j, cons$j_germline_start, cons$j_germline_end)
  if(nchar(cons$germline_alignment) > nchar(uca)){
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
                             as.numeric(cons$np1_length),
                             as.numeric(cons$d_germline_length),
                             as.numeric(cons$np2_length)))){
      diff <- nchar(ref_j) - (nchar(sub$data[[1]]@germline) - sum(nchar(ref_v),  
                                                                  cons$np1_length, cons$d_germline_length, cons$np2_length))
      ref_j <- substring(ref_j, 1, nchar(ref_j) - diff)
    }
  } else{
    if(nchar(ref_j) > (nchar(sub$data[[1]]@lgermline) - sum(nchar(ref_v), 
                             as.numeric(cons$np1_length), 
                             as.numeric(cons$d_germline_length),
                             as.numeric(cons$np2_length)))){
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
      if(!(i %in% v_con_indx)){
        if(length(v_groups[[i]]) == 3){
          codon_value <- paste0(ref_v[v_groups[[i]]], collapse = "")
          if(alakazam::translateDNA(codon_value) != "*"){
            temp <- temp[temp$codon == codon_value,]
          } else{
            # remove the T/U and take the rest -- all stop codons start with T/U with
            ending_values <- substring(codon_value, 2, 3)
            temp <- temp[endsWith(temp$codon, ending_values),]
            warning("A stop codon was detected in the reference for clone ", 
                    cons$clone_id, "which may indicate the improper reference is",
                    " being used")
            writeLines(paste("A stop codon was detected in the reference for clone", 
                              cons$clone_id, "which may indicate the improper",
                             "reference is being used"), 
                       file.path(subDir, "annotation_warning.txt"))
          }
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
      if(!i %in% j_con_indx){
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

# gets the tree table from the tree building process
getTreeTable <- function(build, dir, subDir, clone_ids, repertoire_wide){
  if(build == "igphyml"){
    if(repertoire_wide){
      tree_df <- suppressWarnings(read.table(file = file.path(dir, "sample", "sample_recon_sample",
                                                              paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")), 
                                             header = F, sep = "\t"))
      file.copy(file.path(dir, "sample", "sample_recon_sample",
                          paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")), 
                file.path(subDir, paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")))
    } else{
      tree_df <- suppressWarnings(read.table(file = file.path(subDir, "sample", "sample_recon_sample",
                                                              paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")), 
                                             header = F, sep = "\t"))
      file.copy(file.path(subDir, "sample", "sample_recon_sample",
                          paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")), 
                file.path(subDir, paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")))
    }
  } else if(build == "pml"){
    tree_df <- suppressWarnings(read.table(file = file.path(subDir, "codon_table.txt"), 
                                           header = F, sep = "\t"))
    file.rename(file.path(subDir, "codon_table.txt"), 
                file.path(subDir, paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")))
  }
  
  colnames(tree_df) = c("site", "codon", "partial_likelihood", "log_likelihood_site",
                        "upper_partial_log_likelihood", "upper_partial_likelihood", "equilibrium")
  tree_df$value <- tree_df$partial_likelihood + log(tree_df$equilibrium)
  return(tree_df)
}

# Prepares clones for UCA inference. 
#
# \code{processCloneGermline} Exports two text files that are used as inputs for 
#                             the UCA script.
# @param clone_ids  A clone id for the clones object. This is only used in parallel
# @param clones     A clones object from \link{formatClones}
# @param data       The airr-table associated with the clones object
# @param dir        The directory where data should be saved to 
# @param build      What tree building method to use
# @param id         The run id
# @param resolve_germ Resolve the V and J genes within the clone?
# @param all_germlines The germlines table needed for resolve_germ
# @param v_call     The name of the v annotation column
# @param j_call     The name of the j annotation column 
# @param locus      The name of the locus column 
# @param quiet      How much noise to print out
# @param clone      The name of the proper clone_id column to use
# @param seq_id     The name of the sequence_id column
# @param chain      HL or H?
# @param check_genes Check if the inferred V/J lengths go into the inferred cdr3 region and adjust accordingly. 
# @param references IMGT references read in using \link{readIMGT}
# @param repertoire_wide Were the trees made using repertoire_wide parameters or do they need to made still?
# @param exec       The exec file path for the tree building method (igphyml)
# @param igblast    Exec for igblast
# @param igblast_database path to where the igblast database is 
# @param ref_path   The path to the reference parent folder to use with igblast
# @param igdata     The path to what internal data should be used for igblast
# @param organism   The organism to use igblast with 
# @param partition The partition to use when building the tree
# @param trunklength  The specified trunk length
# @param v_germ_start V germline start column 
# @param v_germ_end   V germline end column 
# @param j_germ_start J germline start column 
# @param j_germ_end   J germline end column 
# @param search       method for the OLGA bit of the UCA process
#
processCloneGermline <- function(clone_ids, clones, data, dir, build, id,
                                 resolve_germ = FALSE, all_germlines = NULL, 
                                 v_call = "v_call", j_call = "j_call", 
                                 locus = "locus", quiet = 0, clone = "clone_id",
                                 seq_id = "sequence_id", chain = "H", 
                                 check_genes = FALSE, references = NULL, 
                                 repertoire_wide = FALSE, exec = NULL, 
                                 igblast = NULL, igdata= NULL,
                                 igblast_database = NULL, ref_path = NULL,
                                 organism = 'human', partition = "single", 
                                 trunklength = NULL, v_germ_start = "v_germline_start", 
                                 v_germ_end = "v_germline_end", 
                                 j_germ_start = "j_germline_start", 
                                 j_germ_end = "j_germline_end", 
                                 search = "codon", ...){
  sub <- dplyr::filter(clones, !!rlang::sym("clone_id") == clone_ids)
  sub_data <- dplyr::filter(data, !!rlang::sym("clone_id") == clone_ids)
  subDir <- file.path(dir, paste0(id, "_",clone_ids))
  if(!dir.exists(subDir)){
    dir.create(subDir, recursive = T)
  }
  if(quiet > 0){
    print("constructing trees")
  }
  if(build == "igphyml"){
    sub <- getTrees(sub, build = build, exec = exec, rm_temp = FALSE, dir = subDir,
                    asrp = TRUE, chain = chain, nproc = 1, partition = partition,
                    trunkl = trunklength, ...)
  } else if(build == "pml"){
    sub <- getTrees(sub, build = build, rm_temp = FALSE, dir = subDir, asrp = TRUE,
                       nproc = 1, ...)
  } else{
    stop("the tree bulding method ", build, "is not supported")
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
  
  tree_df <- getTreeTable(build, dir, subDir, clone_ids, repertoire_wide)
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
  
  if(resolve_germ){
    if(sub$data[[1]]@phylo_seq == "hlsequence"){
      has_multiple <- all_germlines[all_germlines$clone_id == clone_ids,]
      heavy_indx <- grepl("^IGH", has_multiple$v_call)
      saveRDS(has_multiple, file.path(subDir, "all_germlines.rds"))
      has_multiple_light <- has_multiple[!heavy_indx,]
      has_multiple <- has_multiple[heavy_indx,]
      # make sure that has_multiple_light has the padded ungapped 
      has_multiple_light$ungapped <- unlist(lapply(1:nrow(has_multiple_light), function(z){
        value <- has_multiple_light$ungapped[z]
        if(nchar(value) %% 3 != 0){
          padding <- 3 - nchar(value) %% 3
          value <- paste0(value, paste0(rep("N", padding), collapse = ""))
        }
        return(value)
      }))
      saveRDS(has_multiple_light, file.path(subDir, "all_germlines_light.rds"))
      germlines_light <- do.call(rbind, parallel::mclapply(1:nrow(has_multiple_light),
                      function(x){
                         value <- has_multiple_light$ungapped[x]
                         v_alt <- substring(value, 1, nchar(v_light))
                         j_alt <- substring(value, nchar(paste0(v_light, light_cdr3)) + 1, nchar(value))
                         if(nchar(j_alt) > nchar(j_light)){
                           j_alt <- substring(j_alt, 1, nchar(j_light))
                         }
                         sub_tree_df <- tree_df_light[tree_df_light$site %in%
                                                        c(1:(nchar(v_light)/3),
                                                          (max(tree_df_light$new_site) - (nchar(j_light)/3) + 1):
                                                            max(tree_df_light$new_site)),]
                         # get the likelihood of the v_alt and j_alt combo
                         vj <- paste0(v_alt, j_alt, collapse = "")
                         gene_list <- strsplit(vj, "")[[1]]
                         groupedSeq <- split(gene_list, ceiling(seq_along(gene_list) / 3))
                         if("N" %in% groupedSeq[[length(groupedSeq)]]){
                           # add the N option to the J
                           df_row <- sub_tree_df[sub_tree_df$site == max(sub_tree_df$site),]
                           #df_row$value <- df_row$partial_likelihood + log(df_row$equilbrium)
                           value <- sum(df_row$value)
                           new_row <- sub_tree_df[1,]
                           new_row$site <- max(sub_tree_df$site)
                           new_row$codon <- paste0(groupedSeq[[length(groupedSeq)]], collapse = "")
                           new_row$partial_likelihood <- value
                           sub_tree_df <- rbind(sub_tree_df, new_row)
                         }
                         likelihood <- unlist(lapply(1:length(groupedSeq), function(i){
                           codon <- paste0(groupedSeq[[i]], collapse = "")
                           sitedf <- sub_tree_df[sub_tree_df$site == unique(sub_tree_df$site)[i],]
                           return(sitedf$partial_likelihood[sitedf$codon == codon])
                         }))
                         likelihood <- sum(likelihood)
                         temp <- data.frame(clone_id = x, likelihood = likelihood, v = v_alt, j = j_alt,
                                            v_call = has_multiple_light$v_call[x], j_call = has_multiple_light$j_call[x])
                         return(temp)
                       }, mc.cores = nproc))
      index <- which(germlines_light$likelihood == max(germlines_light$likelihood, na.rm = TRUE))[1]
      germlines_light <- germlines_light[index,]
      if("N" %in% strsplit(germlines_light$v, "")[[1]]){
        stop(paste("There is a N found in resolved V gene for clone", clone_ids))
      }
      v_light <- germlines_light$v
      j_light <- germlines_light$j
      saveRDS(germlines_light, file.path(subDir, "most_like_germlines_light.rds"))
    } else{
      v_light <- NULL
      j_light <- NULL
    }
    if(sub$data[[1]]@phylo_seq == "sequence"){
      has_multiple <- all_germlines[all_germlines$clone_id == clone_ids,]
      heavy_indx <- grepl("^IGH", has_multiple$v_call)
      has_multiple <- has_multiple[heavy_indx,]
    } else if(sub$data[[1]]@phylo_seq == "lsequence"){
      has_multiple <- all_germlines[all_germlines$clone_id == clone_ids,]
      light_indx <- !grepl("^IGH", has_multiple$v_call)
      has_multiple <- has_multiple[light_indx,]
    }
    has_multiple$ungapped <- unlist(parallel::mclapply(1:nrow(has_multiple), function(y){
      current <- has_multiple$ungapped[y]
      if(nchar(current) %% 3 > 0){
        current <- paste0(current, paste0(rep("N", 3 - nchar(current) %% 3), collapse = ""))
      }
      return(current)
    }))
    saveRDS(has_multiple, file.path(subDir, "all_germlines.rds"))
    germlines <- do.call(rbind, lapply(1:nrow(has_multiple), function(z){
      missing <- dplyr::setdiff(1:max(sub$data[[1]]@numbers[1:nchar(sub$data[[1]]@germline)]),
                                sub$data[[1]]@numbers[1:nchar(sub$data[[1]]@germline)])
      value <- paste0(strsplit(has_multiple$germline[z], "")[[1]][-missing], collapse = "")
      if(nchar(value) %% 3 != 0){
        value <- paste0(value, paste0(rep("N", (3-nchar(value) %% 3)), collapse = ""))
      }
      v_alt <- substring(value, 1, v_len)
      j_alt <- substring(value, sum(nchar(mrcacdr3), v_len) + 1, 
                         nchar(value))
      j_alt <- substring(j_alt, 1, nchar(j))
      v_stop <- nchar(v)/3
      j_start_new <- nchar(j)/3
      sub_df <- dplyr::filter(tree_df, !!rlang::sym("site") %in% c(0:(v_stop-1)) | 
                                !!rlang::sym("site") %in% c((v_stop + nchar(mrcacdr3)/3): 
                                                              max(tree_df$site)))
      vj <- paste0(v_alt, j_alt, collapse = "")
      gene_list <- strsplit(vj, "")[[1]]
      groupedSeq <- split(gene_list, ceiling(seq_along(gene_list) / 3))
      if("N" %in% groupedSeq[[length(groupedSeq)]]){
        # add the N option to the J
        df_row <- sub_df[sub_df$site == max(sub_df$site),]
        #df_row$value <- df_row$partial_likelihood + log(df_row$equilbrium)
        value <- sum(df_row$value)
        new_row <- sub_df[1,]
        new_row$site <- max(sub_df$site)
        new_row$codon <- paste0(groupedSeq[[length(groupedSeq)]], collapse = "")
        new_row$partial_likelihood <- value
        sub_df <- rbind(sub_df, new_row)
      }
      likelihood <- unlist(lapply(1:length(groupedSeq), function(i){
        codon <- paste0(groupedSeq[[i]], collapse = "")
        sitedf <- sub_df[sub_df$site == unique(sub_df$site)[i],]
        return(sitedf$partial_likelihood[sitedf$codon == codon])
      }))
      likelihood <- sum(likelihood)
      temp <- data.frame(clone_id = z, likelihood = likelihood, v = v_alt, j = j_alt,
                         v_call = has_multiple$v_call[z], j_call = has_multiple$j_call[z])
      return(temp)
    }))
    index <- which(germlines$likelihood == max(germlines$likelihood, na.rm = TRUE))[1]
    germlines <- germlines[index,]
    if("N" %in% strsplit(germlines$v, "")[[1]] || "." %in% strsplit(germlines$v, "")[[1]]){
      stop(paste("There is a N found in resolved V gene for clone", clone_ids))
    }
    v <- germlines$v
    j <- germlines$j
    saveRDS(germlines, file.path(subDir, "most_likely_germlines.rds"))
    if(!is.null(igblast)){
      sub_data <- sub_data[sub_data[[clone]] == sub$clone_id,]
      if(quiet > 0){
        print("updating cuts based on MRCA")
      }
      subdata <- tryCatch({
        updateAIRRGerm(airr_data = sub_data, clones = sub, igblast = igblast, 
                       igblast_database = igblast_database, references = ref_path, 
                       igdata = igdata, organism = organism, locus = "Ig", 
                       outdir = subDir, nproc = 1, clone = clone, v = v, 
                       j = j, v_light = v_light, j_light = j_light, ...)
      }, error = function(e) {
        sub_data
      })
      sub_data <- subdata
    } else{
      if(sub$data[[1]]@phylo_seq == "sequence" || sub$data[[1]]@phylo_seq == "hlsequence"){
        if(nrow(sub_data[sub_data[[clone]] == clone_ids & sub_data[[v_call]] == germlines$v_call
                     & sub_data[[j_call]] == germlines$j_call,]) > 0){
          con <- findConsensus(sub_data[sub_data[[v_call]] == germlines$v_call &
                                          sub_data[[j_call]] == germlines$j_call & 
                                          sub_data[[locus]] == "IGH",], 
                               clone_id = clone_ids)
          con <- sub_data[sub_data[[seq_id]] == con$cons_id,]
          sub_data[[v_call]][sub_data[[locus]] == "IGH"] <- con$v_call
          sub_data[[v_germ_start]][sub_data[[locus]] == "IGH"] <- con$v_germline_start
          sub_data[[v_germ_end]][sub_data[[locus]] == "IGH" ] <- con$v_germline_end
          sub_data[[j_call]][sub_data[[locus]] == "IGH"] <- con$j_call
          sub_data[[j_germ_start]][sub_data[[locus]] == "IGH"] <- con$j_germline_start
          sub_data[[j_germ_end]][sub_data[[locus]] == "IGH"] <- con$j_germline_end
        }
      } else if(sub$data[[1]]@phylo_seq == "lsequence"){
        if(nrow(sub_data[sub_data[[v_call]] == germlines$v_call
                     & sub_data[[j_call]] == germlines$j_call & 
                     sub_data[[locus]] != "IGH",]) > 0){
          con <- findConsensus(sub_data[sub_data[[v_call]] == germlines$v_call &
                                          sub_data[[j_call]] == germlines$j_call,], 
                               clone_id = clone_ids)
          con <- sub_data[sub_data[[seq_id]] == con$cons_id,]
          sub_data[[v_call]][sub_data[[locus]] != "IGH"] <- con$v_call
          sub_data[[v_germ_start]][sub_data[[locus]] != "IGH"] <- con$v_germline_start
          sub_data[[v_germ_end]][sub_data[[locus]] != "IGH"] <- con$v_germline_end
          sub_data[[j_call]][sub_data[[locus]] != "IGH"] <- con$j_call
          sub_data[[j_germ_start]][sub_data[[locus]] != "IGH"] <- con$j_germline_start
          sub_data[[j_germ_end]][sub_data[[locus]] != "IGH"] <- con$j_germline_end
        }
      }
      if(sub$data[[1]]@phylo_seq == "hlsequence"){
        if(nrow(sub_data[sub_data[[v_call]] == germlines$v_call
                     & sub_data[[j_call]] == germlines$j_call & 
                     sub_data[[locus]] != "IGH",]) > 0){
          con <- findConsensus(sub_data[sub_data[[v_call]] == germlines_light$v_call &
                                          sub_data[[j_call]] == germlines_light$j_call,], 
                               clone_id = clone_ids)
          con <- sub_data[sub_data[[seq_id]] == con$cons_id,]
          sub_data[[v_call]][sub_data[[locus]] != "IGH"] <- con$v_call
          sub_data[[v_germ_start]][sub_data[[locus]] != "IGH"] <- con$v_germline_start
          sub_data[[v_germ_end]][sub_data[[locus]] != "IGH"] <- con$v_germline_end
          sub_data[[j_call]][sub_data[[locus]] != "IGH"] <- con$j_call
          sub_data[[j_germ_start]][sub_data[[locus]] != "IGH"] <- con$j_germline_start
          sub_data[[j_germ_end]][sub_data[[locus]] != "IGH"] <- con$j_germline_end
        }
      }
    }
    # update the germline of the clone and rerun the tree building 
    if(quiet > 0){
      print("reconstructing trees")
    }
    
    # update the germline 
    if(sub$data[[1]]@phylo_seq == "sequence"){
      germline_junc <- substring(sub$data[[1]]@germline, nchar(v) + 1,
                                 nchar(v) + nchar(mrcacdr3))
      sub$data[[1]]@germline <- paste0(v, germline_junc, j)
    } else if(sub$data[[1]]@phylo_seq == "lsequence"){
      germline_junc <- substring(sub$data[[1]]@lgermline, nchar(v_light) + 1,
                                 nchar(v_light) + nchar(light_cdr3))
      sub$data[[1]]@lgermline <- paste0(v_light, germline_junc, j_light)
    } else{
      germline_junc <- substring(sub$data[[1]]@germline, nchar(v) + 1,
                                 nchar(v) + nchar(mrcacdr3))
      sub$data[[1]]@germline <- paste0(v, germline_junc, j)
      lgermline_junc <- substring(sub$data[[1]]@lgermline, nchar(v_light) + 1,
                                  nchar(v_light) + nchar(light_cdr3))
      sub$data[[1]]@lgermline <- paste0(v_light, lgermline_junc, j_light)
      sub$data[[1]]@hlgermline <- paste0(v, germline_junc, j, 
                                         v_light, lgermline_junc, j_light)
    }
    # rename the old folder 
    file.rename(file.path(subDir, "sample"), file.path(subDir, "masked_sample"))
    
    if(build == "igphyml"){
      sub <- getTrees(sub, build = build, exec = exec, rm_temp = FALSE, dir = subDir,
                      asrp = TRUE, chain = chain, nproc = 1, partition = partition,
                      trunkl = trunklength, ...)
    } else if(build == "pml"){
      sub <- getTrees(sub, build = build, rm_temp = FALSE, dir = subDir, asrp = TRUE,
                      nproc = 1, ...)
    } else{
      stop("the tree bulding method ", build, "is not supported")
    }
    tree_df <- getTreeTable(build, dir, subDir, clone_ids, repertoire_wide)
  } else{
    if(!is.null(igblast)){
      subdata <- data[data[[clone]] == sub$clone_id,]
      if(quiet > 0){
        print("updating cuts based on MRCA")
      }
      subdata <- tryCatch({
        updateAIRRGerm(airr_data = subdata, clones = sub, igblast = igblast, 
                       igblast_database = igblast_database, references = ref_path, 
                       igdata = igdata, organism = organism, locus = "Ig", 
                       outdir = subDir, nproc = 1, clone = clone, ...)
      }, error = function(e) {
        subdata
      })
      data <- subdata
    }
  }
  saveRDS(sub, file.path(subDir, "clone.rds"))
  if(search %in% c("codon_depth", "nt_depth")){
    if(sub$trees[[1]]$edge_type != "genetic_distance"){
      sub <- scaleBranches(sub, edge_type = "genetic_distance")
    }
    germline_node <- ape::getMRCA(sub$trees[[1]], sub$trees[[1]]$tip.label)
    mrca_node <- ape::getMRCA(sub$trees[[1]], setdiff(sub$trees[[1]]$tip.label, "Germline"))
    co <- ape::dist.nodes(sub$trees[[1]])
    dist <- data.frame(co[germline_node:nrow(co), germline_node])
    colnames(dist) <- "node_dist_grm"
    rownames(dist) <- germline_node:nrow(co)
    indx <- which(rownames(dist) == mrca_node)
    trunk_length <- dist$node_dist_grm[indx]
    if(sub$data[[1]]@phylo_seq != "hlsequence"){
      mrca <- getNodeSeq(sub, mrca_node, sub$trees[[1]])[[1]]
      mrca <- gsub("\\.", "", mrca)
      cdr3_sites <- length(cdr3_index)/3
      center <- trunk_length*cdr3_sites
    } else{
      hmrca <- getNodeSeq(sub, mrca_node, sub$trees[[1]])[[1]]
      lmrca <- getNodeSeq(sub, mrca_node, sub$trees[[1]])[[2]]
      hmrca <- gsub("\\.", "", hmrca)
      lmrca <- gsub("\\.", "", lmrca)
      center <- trunk_length*((length(cdr3_index)-6)/length(heavy_r))
      lcdr3_index <- length(which(light_r == "cdr3"))
      lcdr3_sites <- length(lcdr3_index)/3
      lcenter <- trunk_length*lcdr3_sites
      cdr3_sites <- length(cdr3_index)/3
      center <- trunk_length*cdr3_sites
      file_path_lmrca <- file.path(subDir, paste("lmrca.txt"))
      file_path_lcenter <- file.path(subDir, paste("lcenter.txt"))
      writeLines(paste0(lmrca, collapse = ""), con = file_path_lmrca)
      writeLines(paste(lcenter), con = file_path_lcenter)
    }
    file_path_mrca <- file.path(subDir, paste("mrca.txt"))
    file_path_center <- file.path(subDir, paste("center.txt"))
    writeLines(paste0(mrca, collapse = ""), con = file_path_mrca)
    writeLines(paste(center), con = file_path_center)
  }
  if(check_genes){
    if(sub$data[[1]]@phylo_seq == "hlsequence"){
      heavy_vals <- checkGenesUCA(sub = sub, data = data, v = v, mrcacdr3 = mrcacdr3,
                                  j = j, references = references, tree_df = tree_df,
                                  subDir = subDir, clone_ids = clone_ids, chain = "H",
                                  r = heavy_r, clone = clone)
      light_vals <- checkGenesUCA(sub = sub, data = data, v = v_light, mrcacdr3 = light_cdr3,
                                  j = j_light, references = references, tree_df = tree_df_light,
                                  subDir = subDir, clone_ids = clone_ids, chain = "L",
                                  r = light_r, clone = clone)
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
                                  r = r, clone = clone)
      v <- heavy_vals$v
      mrcacdr3 <- heavy_vals$cdr3
      j <- heavy_vals$j
    } else if(sub$data[[1]]@phylo_seq == "lsequence"){
      light_vals <- checkGenesUCA(sub = sub, data = data, v = v, mrcacdr3 = mrcacdr3,
                                  j = j, references = references, tree_df = tree_df,
                                  subDir = subDir, clone_ids = clone_ids, chain = "L",
                                  r = r, clone = clone)
      v <- light_vals$v
      mrcacdr3 <- light_vals$cdr3
      j <- light_vals$j
    }
  }
  
  if(quiet > 0){
    print(paste("sucessfully obtained most likely junction for", clone_ids))
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
#

callOlga <- function(clones, dir, uca_script, python, max_iters, nproc, id, model_folder,
                     model_folder_igk, model_folder_igl, quiet, search){
  clone_ids <- paste0(unlist(lapply(clones$clone_id, function(z){
    value <- z
    if(clones$data[[which(clones$clone_id == z)]]@phylo_seq == "hlsequence"){
      value <- append(value, z)
    }
    return(value)
  })), collapse = ",")
  starting_germlines <- paste0(unlist(lapply(clones$clone_id, function(z){
    value <- path.expand(file.path(dir, paste0(id, "_", z), "olga_testing_germline.txt"))
    if(clones$data[[which(clones$clone_id == z)]]@phylo_seq == "hlsequence"){
      value <- append(value, path.expand(file.path(dir, paste0(id, "_", z), "olga_testing_germline_light.txt")))
    }
    return(value)
  })), collapse = ",")
  junction_location <- paste0(unlist(lapply(clones$clone_id, function(z){
    value <- path.expand(file.path(dir, paste0(id, "_", z), "olga_junction_positions.txt"))
    if(clones$data[[which(clones$clone_id == z)]]@phylo_seq == "hlsequence"){
      value <- append(value, path.expand(file.path(dir, paste0(id, "_", z), "olga_junction_positions_light.txt")))
    }
    return(value)
  })), collapse = ",")
  tree_tables <- paste0(unlist(lapply(clones$clone_id, function(z){
    if(clones$data[[which(clones$clone_id == z)]]@phylo_seq == "sequence"){
      value <- path.expand(file.path(dir, paste0(id, "_", z), 
                                     paste0(z, ".fasta_igphyml_rootprobs_hlp.txt")))
    } else if(clones$data[[which(clones$clone_id == z)]]@phylo_seq == "hlsequence"){
      value <- path.expand(file.path(dir, paste0(id, "_", z), 
                                     "heavy_table.txt"))
      value <- append(value, path.expand(file.path(dir, paste0(id, "_", z),
                                                   "light_table.txt")))
    } else if(clones$data[[which(clones$clone_id == z)]]@phylo_seq == "lsequence"){
      value <- path.expand(file.path(dir, paste0(id, "_", z), 
                                     paste0(z, ".fasta_igphyml_rootprobs_hlp.txt")))
    }
    return(value)
  })), collapse = ",")
  chains <- paste0(unlist(lapply(clones$clone_id, function(z){
    loci <- strsplit(clones$locus[which(clones$clone_id == z)], ",")[[1]]
  })), collapse = ",")
  if(search %in% c("codon_depth", "nt_depth")){
    mrcas <- paste0(unlist(lapply(clones$clone_id, function(z){
      if(clones$data[[which(clones$clone_id == z)]]@phylo_seq == "sequence"){
        value <- path.expand(file.path(dir, paste0(id, "_", z), "mrca.txt"))
      } else if(clones$data[[which(clones$clone_id == z)]]@phylo_seq == "hlsequence"){
        value <- path.expand(file.path(dir, paste0(id, "_", z), "mrca.txt"))
        value <- path.expand(file.path(dir, paste0(id, "_", z), "lmrca.txt"))
      } else if(clones$data[[which(clones$clone_id == z)]]@phylo_seq == "lsequence"){
        value <- path.expand(file.path(dir, paste0(id, "_", z), "mrca.txt"))
      }
      return(value)
    })), collapse = ",")
    
    centers <- paste0(unlist(lapply(clones$clone_id, function(z){
      if(clones$data[[which(clones$clone_id == z)]]@phylo_seq == "sequence"){
        value <- path.expand(file.path(dir, paste0(id, "_", z), "center.txt"))
      } else if(clones$data[[which(clones$clone_id == z)]]@phylo_seq == "hlsequence"){
        value <- path.expand(file.path(dir, paste0(id, "_", z), "center.txt"))
        value <- path.expand(file.path(dir, paste0(id, "_", z), "lcenter.txt"))
      } else if(clones$data[[which(clones$clone_id == z)]]@phylo_seq == "lsequence"){
        value <- path.expand(file.path(dir, paste0(id, "_", z), "center.txt"))
      }
      return(value)
    })), collapse = ",")
    
  } else{
    mrcas <- "NULL"
    centers <- "NULL"
  }
  
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
    "--search", search, 
    "--mrcas", mrcas, 
    "--centers", centers
  )
  args_keys <- args[seq(1, length(args), 2)]
  args_values <- args[seq(2, length(args), 2)]
  args_list <- stats::setNames(as.list(args_values), gsub("^--", "", args_keys))
  json_str <- paste0(
    "{\n",
    paste(
      sprintf('  "%s": "%s"', names(args_list), gsub('"', '\\"', args_list)),
      collapse = ",\n"
    ),
    "\n}"
  )
  writeLines(json_str, file.path(path.expand(dir), "olga_args.json"))
  cmd <- paste(
    shQuote(python), 
    shQuote(path.expand(uca_script)), 
    "--args_json", shQuote(file.path(path.expand(dir), "olga_args.json"))
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

#' \link{createAllGermlines} Creates all possible germlines for a clone
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
#' @return A data frame with all possible reconstructed germlines
#' @details Return object adds/edits following columns:
#' \itemize{
#'   \item  \code{seq}:  Sequences potentially padded  same length as germline
#'   \item  \code{germline_alignment}: Full length germline
#'   \item  \code{germline_alignment_d_mask}: Full length, D region masked
#'   \item  \code{vonly}:   V gene segment of germline if vonly=TRUE
#'   \item  \code{regions}: String of VDJ segment in position if use_regions=TRUE
#' }
#' @seealso \link{createGermlines} \link{buildAllClonalGermlines}, \link{stitchVDJ}
#' @export
#' 
createAllGermlines <- function(data, references, locus="locus", trim_lengths=FALSE,
                               force_trim=FALSE, nproc=1, seq="sequence_alignment",
                               v_call="v_call", d_call="d_call", j_call="j_call",
                               amino_acid=FALSE, id="sequence_id", clone="clone_id",
                               v_germ_start="v_germline_start", v_germ_end="v_germline_end",
                               v_germ_length="v_germline_length", d_germ_start="d_germline_start",
                               d_germ_end="d_germline_end", d_germ_length="d_germline_length",
                               j_germ_start="j_germline_start", j_germ_end="j_germline_end",
                               j_germ_length="j_germline_length", np1_length="np1_length",
                               np2_length="np2_length", na.rm=TRUE, fields=NULL,
                               verbose=0, ...){
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
  germlines <- do.call(rbind, parallel::mclapply(1:nrow(unique_clones), function(x){
    sub <- dplyr::right_join(data, unique_clones[x,,drop=F], by=c(clone,fields))
    if(verbose > 0){
      print(x)
    }
    sub_germlines <- do.call(rbind, lapply(unique(sub[[locus]]), function(l){
      buildAllClonalGermlines(receptors = sub[sub[[locus]] == l,],
        references = references, chain = l, seq = seq, v_call = v_call, 
        d_call = d_call, j_call = j_call, amino_acid = amino_acid,
        id = id, clone = clone, v_germ_start = v_germ_start,
        v_germ_end = v_germ_end, v_germ_length = v_germ_length,
        d_germ_start = d_germ_start, d_germ_end = d_germ_end, 
        d_germ_length = d_germ_length, j_germ_start = j_germ_start, 
        j_germ_end = j_germ_end, j_germ_length = j_germ_length,
        np1_length = np1_length, np2_length = np2_length, ...)
    }))
  }, mc.cores = nproc))
  return(germlines)
}

#' \link{buildAllClonalGermlines} Determines and builds all possible germlines for a clone
#' @param receptors        AIRR-table containing sequences from one clone
#' @param references       Full list of reference segments, see \link{readIMGT}
#' @param chain            chain in \code{references} being analyzed
#' @param use_regions      Return string of VDJ regions? (optional)
#' @param vonly            Return germline of only v segment?
#' @param seq              Column name for sequence alignment
#' @param id               Column name for sequence ID
#' @param clone            Column name for clone ID
#' @param v_call           Column name for V gene segment gene call
#' @param j_call           Column name for J gene segment gene call
#' @param v_germ_start     Column name of index of V segment start within germline
#' @param v_germ_end       Column name of index of V segment end within germline
#' @param v_germ_length    Column name of index of V segment length within germline
#' @param d_germ_start     Column name of index of D segment start within germline
#' @param d_germ_end       Column name of index of D segment end within germline
#' @param d_germ_length    Column name of index of D segment length within germline
#' @param j_germ_start     Column name of index of J segment start within germline
#' @param j_germ_end       Column name of index of J segment end within germline
#' @param j_germ_length    Column name of index of J segment length within germline
#' @param np1_length       Column name in receptor specifying np1 segment length 
#' @param np2_length       Column name in receptor specifying np2 segment length
#' @param j_germ_aa_length Column name of J segment amino acid length (if amino_acid=TRUE)
#' @param amino_acid       Perform reconstruction on amino acid sequence (experimental)
#' @param threshold        The maximum germline length difference within a clone. Default is 3. 
#' @param ...              Additional arguments passed to \link{buildGermline}
#' @return A data frame with all possible reconstructed germlines
#' @seealso \link{createAllGermlines} \link{buildGermline}, \link{stitchVDJ}
#' @export
#' 
buildAllClonalGermlines <- function(receptors, references,
                                    chain="IGH", use_regions=FALSE, vonly=FALSE,
                                    seq="sequence_alignment", id="sequence_id", clone="clone_id",
                                    v_call="v_call", j_call="j_call",  v_germ_start="v_germline_start",
                                    v_germ_end="v_germline_end", v_germ_length="v_germline_length", 
                                    d_germ_start="d_germline_start", d_germ_end="d_germline_end", 
                                    d_germ_length="d_germline_length", j_germ_start="j_germline_start", 
                                    j_germ_end="j_germline_end", j_germ_length="j_germline_length", 
                                    np1_length="np1_length", np2_length="np2_length",
                                    j_germ_aa_length= "j_germline_aa_length", amino_acid=FALSE, 
                                    threshold = 3, ...){
  
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
  
  v_dict <- unlist(lapply(receptors[[v_call]],function(x)
    alakazam::getAllele(x, strip_d=FALSE, first = FALSE)))
  v_all <- unique(unlist(strsplit(v_dict, ",")))
  j_dict <- unlist(lapply(receptors[[j_call]],function(x)
    alakazam::getAllele(x, strip_d=FALSE, first = FALSE)))
  j_all <- unique(unlist(strsplit(j_dict, ",")))
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
  cons_normal <- receptors[receptors[[id]] == cons_id,]
  
  combinations <- expand.grid(v_all, j_all)
  combo_in_data <- unlist(lapply(1:nrow(combinations), function(x){
    v_loc <- unlist(lapply(1:length(strsplit(receptors[[v_call]], ",")), function(z){
      if(combinations$Var1[x] %in% strsplit(receptors[[v_call]], ",")[[z]]){
        return(z)
      }
    }))
    j_loc <- unlist(lapply(1:length(strsplit(receptors[[j_call]], ",")), function(z){
      if(combinations$Var2[x] %in% strsplit(receptors[[j_call]], ",")[[z]]){
        return(z)
      }
    }))
    common_rows <- intersect(v_loc, j_loc)
    if(length(common_rows) == 0){
      return(FALSE)
    }else{
      return(TRUE)
    }
  }))
  combinations <- combinations[combo_in_data,]
  
  all_germlines <- c()
  for(x in 1:nrow(combinations)){
    v_cons <- as.character(combinations[x, 1])
    j_cons <- as.character(combinations[x, 2])
    v_match <- unlist(lapply(1:length(v_dict), function(z){
      dict_value <- v_dict[z]
      dict_value <- strsplit(dict_value, ",")[[1]]
      return(v_cons %in% dict_value)
    }))
    j_match <- unlist(lapply(1:length(j_dict), function(z){
      dict_value <- j_dict[z]
      dict_value <- strsplit(dict_value, ",")[[1]]
      return(j_cons %in% dict_value)
    }))
    cons_index <- v_match & j_match & seq_len == max_len
    if(sum(cons_index) == 0){
      cons_index <- v_match & j_match
    }
    # return without germline if no consensus can be found 
    if(sum(cons_index) == 0){
      next
    }
    
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
    
    sub_db <- references[[chain]]
    
    positions <- as.numeric(gregexpr("\\.", cons[[seq]])[[1]])
    
    if(length(sub_db) == 0){
      stop(paste("Reference database for",chain,"is empty"))
    }
    
    germlines <- tryCatch(buildGermline(cons, references=sub_db, seq=seq, 
                                        v_call=v_call, j_call=j_call, j_germ_length=j_germ_length,
                                        amino_acid=amino_acid),error=function(e)e)
    if("error" %in% class(germlines)){
      warning(paste("Clone",unique(receptors[[clone]]),"with v and j genes:", v_cons, j_cons),
              "germline reconstruction error.\n",
              germlines)
      temp <- data.frame(clone_id = unique(receptors[[clone]]),
                         clone_id_unique = paste0(unique(receptors[[clone]]), "_", x),
                         v_call = v_cons,
                         d_call = cons$d_call,
                         j_call = j_cons, 
                         v_len = cons$v_germline_length,
                         np1_len = cons$np1_length,
                         d_len = cons$d_germline_length,
                         np2_len = cons$np2_length, 
                         j_len = cons$j_germline_length,
                         max_seq = max_len, 
                         germline = NA,
                         germline_d_mask = NA, 
                         regions = NA, 
                         positions = paste0(positions, collapse = ","))
      all_germlines <- rbind(all_germlines, temp)
      next
    }
    
    temp <- data.frame(clone_id = unique(receptors[[clone]]),
                       clone_id_unique = paste0(unique(receptors[[clone]]), "_", x),
                       v_call = v_cons,
                       d_call = cons$d_call,
                       j_call = j_cons, 
                       v_len = cons$v_germline_length,
                       np1_len = cons$np1_length,
                       d_len = cons$d_germline_length,
                       np2_len = cons$np2_length, 
                       j_len = cons$j_germline_length,
                       max_seq = max_len, 
                       germline = germlines$full,
                       germline_d_mask = germlines$dmask, 
                       regions = germlines$regions, 
                       positions = paste0(positions, collapse = ","))
    all_germlines <- rbind(all_germlines, temp)
  }
  
  if(sum(is.na(all_germlines$germline)) > 0){
    all_germlines <- all_germlines[-which(is.na(all_germlines$germline)),]
  }
  
  all_germlines$ungapped <- unlist(lapply(1:nrow(all_germlines), function(x){
    numbers <- which(strsplit(cons_normal[[seq]], "")[[1]] == ".")
    germ <- strsplit(all_germlines$germline_d_mask[x], "")[[1]]
    germ <- germ[-numbers]
    germ <- paste0(germ, collapse = "")
    return(germ)
  }))
  all_germlines$nchar <- nchar(all_germlines$ungapped)
  if(max(all_germlines$nchar) - min(all_germlines$nchar) >= threshold){
    all_germlines <- all_germlines[-which(all_germlines$nchar < max(all_germlines$nchar - threshold)),]
  }
  return(all_germlines)
}

# \link{maskAmbigousReferenceSites} Determines and builds all possible germlines for a clone
# @param clones AIRR-table that is the output of \link{formatClones}.
# @param all_germlines A data frame with all possible reconstructed germlines.
# @param clone Column name for the clone ID.
# @param nproc Number of cores to use for parallel processing. Default is 1.
# @param chain What locus chain to focus on. H for heavy HL for heavy + light
# @return A data frame with all possible reconstructed germlines.
maskAmbigousReferenceSites <- function(clones, all_germlines, clone = "clone_id", 
                                       nproc = 1, chain = "H"){
  ambig_table <- data.frame(value = c("R", "K", "S", "Y", "M", "W", "B", "H", 
                                      "N", "D", "V"), 
                            combo_l = c(2, 2, 2, 2, 2, 2, 3, 3, 4, 3, 3), 
                            combo = c("A,G", "G,T", "C,G", "C,T", "A,C", 
                                      "A,T", "C,G,T", "A,C,T", "A,C,G,T", 
                                      "A,G,T", "A,C,G"))
  updated_clones <- do.call(rbind, parallel::mclapply(clones[[clone]], function(x){
    sub <- clones[which(clones[[clone]] == x),]
    sub_germs <- all_germlines[all_germlines[[clone]] == x,]
    if(nrow(sub_germs) %in% c(1,0) & chain == "H"| nrow(sub_germs) %in% c(2,0) & chain == "HL"){
      return(sub)
    }
    if(chain == "H"){
      heavy_indx <- grepl("^IGH", sub_germs$v_call)
      sub_germs <- sub_germs[heavy_indx,]
      comp_df <- do.call(cbind, lapply(1:nrow(sub_germs), function(z){
        temp <- data.frame(strsplit(sub_germs$ungapped[z], "")[[1]])
      }))
      colnames(comp_df) <- sub_germs$clone_id_unique
    } else if(chain == "L"){
      light_indx <- !grepl("^IGH", sub_germs$v_call)
      sub_germs <- sub_germs[light_indx,]
      comp_df <- do.call(cbind, lapply(1:nrow(sub_germs), function(z){
        temp <- data.frame(strsplit(sub_germs$ungapped[z], "")[[1]])
      }))
      colnames(comp_df) <- sub_germs$clone_id_unique
    } else if(chain == "HL"){
      # get the concatenated germline
      heavy_indx <- grepl("^IGH", sub_germs$v_call)
      heavy <- sub_germs[heavy_indx,]
      light <- sub_germs[!heavy_indx,]
      new_germs <- do.call(rbind, lapply(1:nrow(heavy), function(z){
        sub_heavy <- heavy[z,]
        if(nchar(sub_heavy$ungapped) %% 3 != 0){
          padding <- 3 - nchar(sub_heavy$ungapped) %% 3
          sub_heavy$ungapped <- paste0(sub_heavy$ungapped, paste0(rep("N", padding), 
                                                                  collapse = ""))
        }
        if(nrow(light) > 0){
          new_germ_option <- do.call(rbind, lapply(1:nrow(light), function(y){
            temp <- data.frame(clone_id = sub_heavy$clone_id,
                               clone_id_unique = paste0(sub_heavy$clone_id_unique, ",", light$clone_id_unique[y]),
                               v_call = paste0(sub_heavy$v_call, ",", light$v_call[y]),
                               j_call = paste0(sub_heavy$j_call, ",", light$j_call[y]),
                               max_seq = paste0(sub_heavy$max_seq, ",", light$max_seq[y]),
                               germline = paste0(sub_heavy$germline, light$germline[y]), 
                               germline_d_mask = paste0(sub_heavy$germline_d_mask,light$germline_d_mask[y]),
                               regions  = paste0(sub_heavy$regions, ",", light$regions[y]),
                               positions = paste0(sub_heavy$positions,",", light$positions[y]),
                               ungapped = paste0(sub_heavy$ungapped, light$ungapped[y]),
                               nchar = sum(sub_heavy$nchar, light$nchar[y]))
          }))
        } else{
          new_germ_option <- data.frame(clone_id = sub_heavy$clone_id,
                                        clone_id_unique = paste0(sub_heavy$clone_id_unique, ",germline"),
                                        v_call = paste0(sub_heavy$v_call, ",germline"),
                                        j_call = paste0(sub_heavy$j_call, ",germline"),
                                        max_seq = paste0(sub_heavy$max_seq, ",germline"),
                                        germline = paste0(sub_heavy$germline, ",germline"), 
                                        germline_d_mask = paste0(sub_heavy$germline_d_mask,",germline"),
                                        regions  = paste0(sub_heavy$regions, ",germline"),
                                        positions = paste0(sub_heavy$positions,",germline"),
                                        ungapped = paste0(sub_heavy$ungapped, sub$data[[1]]@lgermline),
                                        nchar = sum(sub_heavy$nchar, nchar(sub$data[[1]]@lgermline)))
        }
        
        return(new_germ_option)
      }))
      # then do the comp df
      comp_df <- do.call(cbind, lapply(1:nrow(new_germs), function(z){
        temp <- data.frame(strsplit(new_germs$ungapped[z], "")[[1]])
      }))
      colnames(comp_df) <- new_germs$clone_id_unique
    }
    comp_df$diff <- apply(comp_df, 1, function(row) length(unique(row)) > 1)
    if(nrow(comp_df) %% 3 > 0){
      padding <- 3 - (nrow(comp_df) %% 3)
      for(i in 1:padding){
        temp <- comp_df[nrow(comp_df),]
        temp[1,] <- "N"
        temp$diff <- FALSE
        comp_df <- rbind(comp_df, temp)
      }
    }
    
    if(chain == "H"){
      clone_germ <- strsplit(sub$data[[1]]@germline, "")[[1]]
    } else if(chain == "HL"){
      clone_germ <- strsplit(sub$data[[1]]@hlgermline, "")[[1]]
    } else{
      clone_germ <- strsplit(sub$data[[1]]@lgermline, "")[[1]]
    }
    masked_sites <- which(comp_df$diff)
    for(i in masked_sites){
      row_vals <- comp_df[i,]
      row_vals <- row_vals[,-which(colnames(row_vals) == "diff")]
      unique_vals <- unique(unlist(row_vals))
      combo_string <- paste(sort(unique_vals), collapse = ",")
      result <- ambig_table$value[ambig_table$combo == combo_string]
      if(length(result) == 0){
        result <- "N"
      }
      clone_germ[i]  <- result
    }
    
    if(chain == "H"){
      sub$data[[1]]@germline <- paste0(clone_germ, collapse = "")
    } else if(chain == "HL"){
      sub$data[[1]]@hlgermline <- paste0(clone_germ, collapse = "")
    } else{
      sub$data[[1]]@lgermline <- paste0(clone_germ, collapse = "")
    }
    return(sub)
  }, mc.cores = nproc))
  return(updated_clones)
}

#' \link{getTreesAndUCAs} Construct trees and infer the UCA
#' @param clones        AIRR-table containing sequences \link{formatClones}
#' @param data          The AIRR-table that was used to make the clones object.
#' @param dir           The file path of the directory of where data is saved. NULL is default.
#' @param build         Name of the tree building method. Currently only IgPhyML is supported.
#' @param exec          File path to the tree building executable
#' @param resolve_germ  Resolve the V and J gene annotations within each clone?
#' @param repertoire_wide Build build trees using parameters inferred from the entire dataset?
#' @param partition     The partition model to use with IgPhyML. "single" is the default.
#' @param model_folder  The file path to the OLGA default model files for heavy chains
#' @param model_folder_igk  The file path to the OLGA default model files for IGK
#' @param model_folder_igl  The file path to the OLGA default model files for IGL
#' @param python        Specify the python call for your system. This is the call
#'                      on command line that issues the python you want to use. 
#'                      "python3" by default. 
#' @param id            The run ID
#' @param max_iters     The maximum number of iterations to run before ending
#' @param nproc         The number of cores to use 
#' @param rm_temp       Remove the generated files?
#' @param chain         Set to HL to use both heavy and light chain sequences
#' @param quiet         Amount of noise to print out
#' @param references    Reference genes. See \link{readIMGT}
#' @param clone         The name of the clone id column used in \link{formatClones}
#' @param cell          The name of the cell id in the AIRR table used to generate \link{formatClones}
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
#' @param igdata        The file path to the data used for igblast (NULL for internal data with igblast)
#' @param organism      The type of organism to align to if using igblast. 
#' @param trunklength   The specified trunklength
#' @param ...           Additional arguments passed to various other functions like \link{getTrees} and \link{buildGermline}
#' @return An \code{airrClone} object with trees and the inferred UCA
#' @details Return object adds/edits following columns:
#' \itemize{
#'   \item  \code{trees}:  The phylogenies associated with each clone
#'   \item  \code{UCA}:    The inferred UCA
#' }
#' @seealso \link{getTrees} 
#' @export
getTreesAndUCAs <- function(clones, data, dir = NULL, build = "igphyml", 
                            exec = NULL, resolve_germ = FALSE, 
                            repertoire_wide = FALSE, partition = "single", 
                            model_folder, model_folder_igk = NULL, 
                            model_folder_igl = NULL, python = "python3", 
                            id = "sample", max_iters = 100, nproc = 1, 
                            rm_temp = TRUE, quiet = 0, chain = "H",
                            references = NULL, clone = "clone_id", 
                            cell = "cell_id", sampling_method = 'random', 
                            subsample_size = NA, search = "codon", 
                            check_genes = TRUE, igblast = NULL, 
                            igblast_database = NULL, ref_path = NULL, 
                            igdata = NULL, organism = 'human', 
                            trunklength = NULL, ...){
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
    stop("The python executable provided, ", python, " cannot be executed.")
  }
  
  uca_script <- system.file("get_UCA.py", package = "dowser")
  if (!file.exists(path.expand(uca_script))) {
    stop(
      paste(
        "The required Python script 'get_UCA.py' could not be found or accessed at:",
        path.expand(uca_script), "\n",
        "This script should be included with the 'dowser' R package.",
        "Please update or reinstall the 'dowser' package to ensure all scripts are present."
      )
    )
  }
  
  if(check_genes & is.null(references)){
    stop("check_genes cannot be run without references. Pass a references object using references call.",
         "References need to be read in using dowser::readIMGT()")
  }
  
  if(resolve_germ & is.null(data) | resolve_germ & is.null(references)){
    stop('resolve_germ requires the data object and references', 
         "References need to be read in using dowser::readIMGT()")
  }
  
  if(!is.null(igblast) & any(sapply(list(igblast_database, ref_path), is.null))){
    stop('if you are planning on updating the germline references using igblast, the following varaibles cannot be NULL:',
         "\n",
         'igblast: the file path to igblast', "\n",
         'igblast_database: the file path to the databases set up for igblast',"\n",
         'ref_path: the path to the imgt parent folder')
  }
  
  if(build == "igphyml"){
    bad_germs <- unlist(parallel::mclapply(1:nrow(clones), function(x){
      germ <- clones$data[[x]]@germline
      if("*" %in% strsplit(translateDNA(germ), "")[[1]]){
        return(clones$clone_id[x])
      }
    }, mc.cores = nproc))
    if(length(bad_germs) > 0){
      warning('Stop codon was detected in the germline of clones ', bad_germs, 
              ' and will be removed from this analysis')
    }
    clones <- clones[!clones$clone_id %in% bad_germs,]
  }
  
  if(resolve_germ){
    all_germlines <- createAllGermlines(data = data, references = references,
                                        nproc = nproc, clone = clone, 
                                        trim_lengths = TRUE, ...)
    saveRDS(all_germlines, file.path(dir, "all_germlines.rds"))
    clones <- maskAmbigousReferenceSites(clones = clones, all_germlines = all_germlines,
                                         nproc = nproc, clone = clone, chain = chain)
  } else{
    all_germlines <- NULL
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
    if(subsample_size == 1 & !is.na(subsample_size)){
      cells <- unlist(lapply(clones$data, function(x) x@data$sequence_id))
      if (!is.na(cell) & cell %in% colnames(data)) {
        filtered <- data[data$sequence_id %in% cells, ]
        cells <- data$sequence_id[data[[cell]] %in% filtered[[cell]]]
      }
      sub_data <- data[data$sequence_id %in% cells,]
      clones <- formatClones(sub_data, nproc = nproc, 
                             filterstop = TRUE, chain = chain, minseq = 1, 
                             dup_singles = T, ...)
    }
    saveRDS(clones, file = file.path(dir, "clones.rds"))
  }
  if(repertoire_wide){
    if(quiet > 0){
      print("constructing trees")
    }
    if(build == "igphyml"){
      if(chain == "HL" & partition != "hl"){
        warning("Paired analysis is being requested but the paired partition is not being requested. To build the best paired trees use partition = 'hl'")
      } 
      clones <- getTrees(clones, build = build, exec = exec, rm_temp = FALSE, dir = dir,
                         asrp = TRUE, chain = chain, nproc = nproc, partition = partition,
                         trunkl = trunklength, ...)


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
      if(is.null(igdata)){
        igdata <- igblast
      }
      data <- updateAIRRGerm(airr_data = data, clones = clones, igblast = igblast, 
                             igblast_database = igblast_database, references = ref_path, 
                             igdata = igdata, organism = organism, locus = 'Ig',
                             outdir = dir, nproc = nproc, clone = clone, ...)
    }
  }
  
  if(quiet > 0){
    print("preparing the clones for UCA analysis")
  }
  clones <- invisible(do.call(rbind, parallel::mclapply(clones$clone_id, function(x){
    processCloneGermline(clone_ids = x, clones = clones, data = data, dir = dir,
                         build = build, id = id, resolve_germ = resolve_germ, 
                         all_germlines = all_germlines, quiet = quiet, 
                         clone = clone, chain = chain, check_genes = check_genes, 
                         exec = exec, references = references, 
                         repertoire_wide = repertoire_wide, igblast = igblast, 
                         igblast_database = igblast_database, ref_path = ref_path, 
                         organism = organism, partition = partition,
                         trunklength = trunklength, search = search, ...) 
    # the problem is happening when rerunning igblast for the MRCA -- see why
  }, mc.cores = nproc)))
  saveRDS(clones, file.path(dir, "clones.rds"))
  # run the UCA
  if(quiet > 0){
    print("running UCA analysis")
  }
  callOlga(clones = clones, dir = dir, model_folder = model_folder,
           model_folder_igk = model_folder_igk, model_folder_igl = model_folder_igl,
           uca_script = uca_script, python = python, max_iters = max_iters,
           nproc = nproc, id = id, quiet = quiet, search = search, ...)

  if(quiet > 0){
    print("updating clones")
  }
  clones <- updateClone(clones, dir, id, nproc)
  saveRDS(clones, file.path(dir, "clones.rds"))
  unlink(rm_dir,recursive=TRUE)
  return(clones)
}
