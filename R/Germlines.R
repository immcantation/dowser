## Functions for constructing clonal germline sequences
## Based closely on CreateGermlines.py
#' \code{downloadIMGT} Download IMGT GENE-DB databases
#' 
#' Loads all reference germlines from an Immcantation-formatted IMGT database.
#' 
# TODO: make auto-download or internal IMGT database
#' @param dir      directory to contain database
#' @param organism  vector of species to download: human, mouse, rat, rabbit, rhesus_monkey
#' @param timeout  max time allowed for each download (temporarily sets global timeout option)
#' @return directory of downloaded IMGT BCR and TCR reference files
#' @details 
#' Recoding of this script: 
#' https://bitbucket.org/kleinstein/immcantation/src/master/scripts/fetch_imgtdb.sh
#' @export
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

#' Write a fasta file of sequences
#' \code{readFasta} reads a fasta file
#' @param    df        dataframe of sequences
#' @param    id        Column name of sequence ids
#' @param    seq       Column name of sequences
#' @param    file      FASTA file for output
#' @param    imgt_gaps Keep IMGT gaps if present?
#' @param    columns   vector of column names to append to sequence id
#'
#' @return   File of FASTA formatted sequences
#' @export
dfToFasta <- function(df, file, id="sequence_id", seq="sequence",
  imgt_gaps=FALSE, columns=NULL){
  if(!"data.frame" %in% class(df)){
    stop("df must be a data.frame or tibble")
  }
  if(!id %in% names(df)){
    stop(id, " column not found in df")
  }
  if(!seq %in% names(df)){
    stop(seq, " column not found in df")
  }

  if(!imgt_gaps){
    seqs <- gsub("\\.","",df[[seq]])
  }else{
    seqs <- df[[seq]]
  }

  if(is.null(columns)){
    ids <- df[[id]]
  }else{
    if(sum(!columns %in% names(df)) > 0){
      nf <- columns[!columns %in% names(df)]
      stop(paste(nf, collapse=","), " not found in df")
    }
    ids <- df[[id]]
    for(column in columns){
      values <- paste0("|", column, "=", df[[column]])
      ids <- paste0(ids, values)
    }
  }
  lines <- paste0(">", ids, "\n", seqs)
  writeLines(lines, con=file)
}

#' \code{imgtToIgblast} Format IMGT database to Igblast database
#' 
#' Cleans IMGT database and creates blast dbs for IgBlast
#' 
#' @param imgt  directory containing IMGT sequences ("dir" in downloadIMGT)
#' @param outdir   output directory for IgBlast databases, will contain /fasta and /database
#' @param igblast  IgBlast home directory location
#' @param organism   organism species to use (must be human, mouse, or rhesus_monkey)
#' @return outdir
#' @details 
#' Recoding of this script: 
#' https://bitbucket.org/kleinstein/immcantation/src/master/scripts/imgt2igblast.sh
#' @export
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

#' \code{assignGenes} Runs IgBlast on a fasta file
#'  
#' @param file      Fasta file of Ig or TR sequences
#' @param igblast   Location of IgBlast program directory, containing /bin/igblastn
#' @param refs     Reference directory of sequences (see downloadIMGT)
#' @param igdata  Internal data directory for IgBlast (defaults to \code{igblast}/internal_data)
#' @param organism    Organism (human, mouse, rhesus_monkey)
#' @param domain_system Currently only IMGT supported
#' @param outfile       Name of AIRR rearrangement file (must end in TSV)
#' @param nproc     Nummber of threads to use
#' @param db_prefix   File prefix for reference fastas
#' @param locus     Ig or TR
#' @param set_igdata  Set IGDATA environment variable?
#' @param return    Return data.frame of output?
#' @param verbose   Print extra info?
#' @return AIRR-rearrangement formatted data frame
#' @details 
#' Runs IgBlast, similar to AssignGenes.py in Changeo.
#' Must have IgBlast downloaded, precompiled binaries recommended:
#' https://ncbi.github.io/igblast/
#' 
#' Note for M1/M2 Mac OS, may be necessary to install IgBlast .dmg and 
#' manually set igdata. Otherwise most issues stem from the internal_data
#' directory location.
#' 
#' @export
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

  if(!is.null(igdata) && set_igdata){
    if(is.null(igdata)){
      print("igdata not specified, using igblast for IGDATA")
      igdata <- igblast
    }

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

#' \code{addGaps} Add IMGT gaps to IgBlast output
#'  
#' @param db      AIRR-rearrangment formatted dataframe outputted from \link{assignGenes}
#' @param gapdb   Root directory of gapped fasta files (\code{dir} in \link{downloadIMGT})
#' @param organism  Organism (human, mouse, rhesus_monkey)
#' @param locus   Ig or TR
#' @param gapped_d  Include IMGT gaps in D regions? Only applicable to 
#' @return AIRR-rearrangement formatted data frame in which sequence_alignment and germline_alignment
#' have been updated with IMGT gaps
#' @details 
#' Similar functionality to MakeDb.py in Change-O. 
#' 
#' @export
addGaps <- function(db, gapdb, organism="human", locus="Ig", gapped_d=FALSE){

  if(!organism %in% c("human", "mouse", "rhesus_monkey")){
    stop(paste("organism must be either human, mouse, rhesus_monkey"))
  }

  if(!locus %in% c("Ig", "TR")){
    stop(paste("locus must be either Ig or TR"))
  }

  gap_files <- list.files(file.path(gapdb, organism, "vdj"))
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
    results <- bind_rows(results, row)
  }
  return(results)
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
#' segment refers to gene segment caegory (V, D, or J)
#' @details Input directory must be formatted to Immcantation standard.
#' See https://changeo.readthedocs.io/en/stable/examples/igblast.html for example
#' of how to download.
#' @examples
#' # vdj_dir contains a minimal example of reference germlines 
#' # (IGHV3-11*05, IGHD3-10*01 and IGHJ5*02)
#' # which are the gene assignments for ExamapleDb[1,]
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
#' instead of nulecotides
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
                    "germline reconstruction error."))
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
#' @param trim_n        Remove trailing Ns from \code{seq} column?
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
createGermlines <- function(data, references, locus="locus", trim_n=FALSE,
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
  if(sum(is.na(data[[clone]])) > 0){
    stop("NA values in clone id column found, please remove.")
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
  trailing_ns <- grepl("N+$", data[[seq]])
  if(sum(trailing_ns) > 0){
    trailing_clones <- dplyr::n_distinct(data[trailing_ns,]$clone_id)
    if(!trim_n){
      warning("Trailing Ns in ", sum(trailing_ns)," sequences in ",
       trailing_clones, 
       " clones, consider setting trim_n=TRUE if germline reconstruction fails")
    }else{
      data[[seq]] <- gsub("N+$", "", data[[seq]])
      cat("Removed trailing Ns from", sum(trailing_ns),"sequences in",
       trailing_clones, "clones\n")
    }
  }

  unique_clones <- unique(data[,unique(c(clone,fields)),drop=F])
  data[['tmp_row_id']] <- 1:nrow(data)
  complete <- parallel::mclapply(1:nrow(unique_clones), function(x){
    sub <- right_join(data, unique_clones[x,,drop=F], by=c(clone,fields))
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
