## Functions for constructing clonal germline sequences
# Based closely on CreateGermlines.py
# Download IMGT GENE-DB databases
# 
# Loads all reference germlines from an Immcantation-formatted IMGT database.
# 
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

    # remove weird list error that sometimes comes up
    glines <- lapply(glines, function(x){
      lcols <- names(which(sapply(names(x), function(y)class(x[[y]])) == "list"))
      for(lcol in lcols){
        x[[lcol]] <- unlist(x[[lcol]])
      }
      x
    })
    glined <- tryCatch(dplyr::bind_rows(glines),
      error=function(e){
        saveRDS(glines, "glines_error.rds")
        stop("caught bind_rows(glines) error")
      })
    glined
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
# @param receptors        An airr table for a specific clone
# @param v_call           The column name for v_call in receptors
# @param j_call           The column name for j_call in receptors
# @param seq              The column name for sequence_alignment in receptors
# @param id               The column name for sequence_id in receptors
# @param v_germ_length    The column name for v_germ_length in receptors
# @param d_germ_length    The column name for d_germ_length in receptors
# @param j_germ_length    The column name for j_germ_length in receptors
# @param np1_length       The column name for np1_length in receptors
# @param np2_length       The column name for np2_length in receptors
findConsensus <- function(receptors, v_call = "v_call", j_call = "j_call",
                          seq = "sequence_alignment", id = "sequence_id", 
                          v_germ_length="v_germline_length", 
                          d_germ_length="d_germline_length", 
                          j_germ_length="j_germline_length", 
                          np1_length="np1_length", np2_length="np2_length"){
  
  if (nrow(receptors) == 0) {
    stop("findConsensus received an empty data frame. ",
         "No sequences were found for this locus/clone combination. ",
         "Check that locus column values and clone membership are correct.")
  }
  
  v_dict <- c()
  j_dict <- c()
  pad_char <- "N"
  
  receptors[[seq]] <- unlist(lapply(seq_len(nrow(receptors)), function(x) {
    if (is.na(receptors[[np1_length]][x]))   receptors[[np1_length]][x]   <- 0
    if (is.na(receptors[[d_germ_length]][x])) receptors[[d_germ_length]][x] <- 0
    if (is.na(receptors[[np2_length]][x]))   receptors[[np2_length]][x]   <- 0
    
    germline_len <- sum(receptors[[v_germ_length]][x],
                        receptors[[np1_length]][x],
                        receptors[[d_germ_length]][x],
                        receptors[[np2_length]][x],
                        receptors[[j_germ_length]][x],
                        na.rm = TRUE) 
    
    seq_len <- nchar(receptors[[seq]][x])
    
    if (is.na(germline_len) || is.na(seq_len)) {
      warning("NA detected in germline or sequence length for row ", x,
              " of clone ", receptors[id][x],
              ". Returning sequence as-is.")
      return(receptors[[seq]][x])
    }
    
    if (germline_len != seq_len) {
      diff  <- abs(germline_len - seq_len)
      value <- substring(receptors[[seq]][x], 1, nchar(receptors[[seq]][x]) - diff)
    } else {
      value <- receptors[[seq]][x]
    }
    return(value)
  }))
  
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
  temp <- data.frame(clone_id = receptors$clone_id[1], locus = receptors$locus[1], v_call = rec_v, 
                     j_call = rec_j, cons_id = cons_id)
  return(temp)
}


# checkGenesUCA is what is run if check_genes = TRUE in getTreesAndUCAs
# @param sub    A clones object for only 1 clone
# @param cons   The airr table entry for the consensus sequence within a clone 
# @param v      A string for the v gene
# @param cdr3   A string of the junction (conserved site to conserved site)
# @param j      A string for the j gene 
# @param tree_df  The df asscoiated with the partial likelihood in tree building 
# @param subDir   The specific directory for the clone of interest 
# @param clone_ids  The specific clone_id of interest 
# @param regions    The regions of sub -- if sub is run with chain == "HL" this 
#                   vector should only have one chain's worth of regions
# @param chain    The chain option used in formatClones ("H", "HL", or "L")
# @param clone    The name of the clone_id varaible in the airr table used to 
#                 create the clone
# @param seq_id   The name of the sequence_id varaible in the airr table used to 
#                 create the clone
checkGenesUCA <- function(sub, cons, v, cdr3, j, tree_df, subDir, clone_ids, regions, 
                          chain = "H", clone = "clone_id", 
                          seq_id = "sequence_id"){
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
  
  if(chain == "H"){
    germline <- sub$data[[1]]@germline
  } else{
    germline <- sub$data[[1]]@lgermline
  }
  
  uca <- strsplit(germline, "")[[1]]
  for(pos in gaps) {
    uca <- append(uca, ".", after = pos - 1)
  }
  uca <- paste(uca, collapse = "")
  for(pos in gaps) {
    regions <- append(regions, "gap", after = pos - 1)
  }
  
  ref_v <- substring(uca, cons$v_germline_start, cons$v_germline_end)
  ref_v <- sub("N.*$", "", ref_v)
  
  if(chain == "L" | sub$data[[1]]@phylo_seq == "lsequence"){
    pad_length <- sum(strsplit(substring(sub$data[[1]]@lgermline, (nchar(sub$data[[1]]@lgermline)-2),
                                         (nchar(sub$data[[1]]@lgermline))), "")[[1]] == "N")
  } else{
    pad_length <- sum(strsplit(substring(sub$data[[1]]@germline, (nchar(sub$data[[1]]@germline)-2),
                                         (nchar(sub$data[[1]]@germline))), "")[[1]] == "N")
  }
  
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
  
  ref_j <- substring( uca, nchar(uca) - cons$j_germline_length + 1 - pad_length,
                      nchar(uca))
  ref_j <- substring(ref_j, 1, nchar(ref_j) - pad_length)
  ref_j <- sub("^.*N", "", ref_j)
  ref_j <- paste0(ref_j, paste(rep("N", pad_length), collapse = ""))

  ref_v <- paste0(strsplit(ref_v, "")[[1]][-gaps], collapse = "")
  j_gaps <- gaps[gaps >= nchar(uca) - nchar(ref_j) + 1]
  
  if(length(j_gaps) > 0){
    j_gaps <- j_gaps - (nchar(uca) - nchar(ref_j) + 1)
    ref_j <- paste0(strsplit(ref_j, "")[[1]][-j_gaps], collapse = "")
  }
  
  if(cons$locus == "IGH"){
    if(nchar(ref_j) > (nchar(sub$data[[1]]@germline) - sum(
      nchar(ref_v), as.numeric(cons$np1_length), as.numeric(cons$d_germline_length),
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
  
  # find the conserved sites to ensure those get saved no matter what
  # V conserved
  v_cons <- (nchar(v) + 1):(nchar(v) + 3)
  # J conserved 
  j_con <- (nchar(v) + nchar(cdr3) - 2):(nchar(v) + nchar(cdr3))
  
  v_con_in_v_indx <- which(sapply(v_groups, function(x) any(x %in% v_cons)))
  j_con_in_v_indx <- which(sapply(v_groups, function(x) any(x %in% j_con)))
  
  v_df <- do.call(rbind, lapply(1:length(v_groups), function(i){
    temp <- tree_df[tree_df$site == i - 1,]
    
    is_v_con <- i %in% v_con_in_v_indx   # V conserved C site
    is_j_con <- i %in% j_con_in_v_indx   # J conserved W/F site (rare edge case)
    
    if(is_v_con){
      # Allow reference codon OR any codon that translates to C
      if(length(v_groups[[i]]) == 3){
        temp <- temp[
          (temp$codon == paste0(ref_v[v_groups[[i]]], collapse = "")) |
            (alakazam::translateDNA(temp$codon) == "C"), ]
        if(alakazam::translateDNA(paste0(ref_v[v_groups[[i]]], collapse = "")) == "C"){
          temp <- temp[temp$codon == paste0(ref_v[v_groups[[i]]], collapse = ""), ]
        }
      } else{
        values <- paste0(ref_v[v_groups[[i]]], collapse = "")
        temp <- temp[(startsWith(temp$codon, values)) |
                       (alakazam::translateDNA(temp$codon) == "C"), ]
      }
    } else if(is_j_con){
      # J conserved site landed in V region: allow reference codon OR W/F
      if(length(v_groups[[i]]) == 3){
        temp <- temp[
          (temp$codon == paste0(ref_v[v_groups[[i]]], collapse = "")) |
            (alakazam::translateDNA(temp$codon) %in% c("W", "F")), ]
        if(alakazam::translateDNA(paste0(ref_v[v_groups[[i]]], collapse = "")) %in% c("W", "F")){
          temp <- temp[temp$codon == paste0(ref_v[v_groups[[i]]], collapse = ""), ]
        }
      } else{
        values <- paste0(ref_v[v_groups[[i]]], collapse = "")
        temp <- temp[(startsWith(temp$codon, values)) |
                       (alakazam::translateDNA(temp$codon) %in% c("W", "F")), ]
      }
    } else{
      # Non-conserved V site, but a conserved site exists nearby: N/stop-aware filter
      if(length(v_groups[[i]]) == 3){
        codon_value <- paste0(ref_v[v_groups[[i]]], collapse = "")
        if(!alakazam::translateDNA(codon_value) %in% c("*", "X")){
          temp <- temp[temp$codon == codon_value, ]
        } else{
          if("N" %in% ref_v[v_groups[[i]]]){
            pattern <- ref_v[v_groups[[i]]]
            non_n_positions <- which(pattern != "N")
            condition <- rep(TRUE, nrow(temp))
            for (pos in non_n_positions) {
              condition <- condition & (substr(temp$codon, pos, pos) == pattern[pos])
            }
            matching_codons <- temp[condition, ]
            value <- sum(matching_codons$value)
            new_row <- temp[1, ]
            new_row$codon <- paste0(ref_v[v_groups[[i]]], collapse = "")
            new_row$partial_likelihood <- value
            temp <- new_row
          } else{
            ending_values <- substring(codon_value, 2, 3)
            temp <- temp[endsWith(temp$codon, ending_values), ]
            warning("A stop codon was detected in the reference for clone ", 
                    cons$clone_id, "which may indicate the improper reference is",
                    " being used")
            writeLines(paste("A stop codon was detected in the reference for clone", 
                             cons$clone_id, "which may indicate the improper",
                             "reference is being used"), 
                       file.path(subDir, "v_annotation_warning.txt"))
          }
        }
      } else{
        if(!"N" %in% ref_v[v_groups[[i]]]){
          values <- paste0(ref_v[v_groups[[i]]], collapse = "")
          temp <- temp[startsWith(temp$codon, values), ]
        } else{
          pattern <- ref_v[v_groups[[i]]]
          non_n_positions <- which(pattern != "N")
          condition <- rep(TRUE, nrow(temp))
          for (pos in non_n_positions) {
            condition <- condition & (substr(temp$codon, pos, pos) == pattern[pos])
          }
          matching_codons <- temp[condition, ]
          value <- sum(matching_codons$value)
          new_row <- temp[1, ]
          new_row$codon <- paste0(ref_v[v_groups[[i]]], collapse = "")
          new_row$partial_likelihood <- value
          temp <- new_row
        }
      }
    }
    return(temp)
  }))
  
  j_df <- tree_df[tree_df$site %in% tail(sort(unique(tree_df$site)),
                                         length(j_groups)), ]
  j_df$new_site <- j_df$site - (min(j_df$site) - 1)
  
  if(chain == "L" | sub$data[[1]]@phylo_seq == "lsequence"){
    offset <- nchar(sub$data[[1]]@lgermline) - max(j_groups[[length(j_groups)]])
  } else{
    offset <- nchar(sub$data[[1]]@germline) - max(j_groups[[length(j_groups)]])
  }
  
  j_groups_num <- lapply(j_groups, function(x) x + offset)
  
  v_con_in_j_indx <- which(sapply(j_groups_num, function(x) any(x %in% v_cons)))
  j_con_in_j_indx <- which(sapply(j_groups_num, function(x) any(x %in% j_con)))
  
  j_df <- do.call(rbind, lapply(1:length(j_groups), function(i){
    temp <- j_df[j_df$new_site == i, ]
    
    is_j_con <- i %in% j_con_in_j_indx   
    is_v_con <- i %in% v_con_in_j_indx  
    
    if(is_j_con){
      # Allow reference codon OR any codon translating to W/F
      if(length(j_groups[[i]]) == 3){
        if(sum("N" %in% ref_j[j_groups[[i]]]) == 0){
          if(alakazam::translateDNA(paste0(ref_j[j_groups[[i]]], collapse = "")) != "*"){
            temp <- temp[temp$codon == paste0(ref_j[j_groups[[i]]], collapse = "") |
                           (alakazam::translateDNA(temp$codon) %in% c("F", "W")), ]
          } else{
            temp <- temp[alakazam::translateDNA(temp$codon) %in% c("F", "W"), ]
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
                       (alakazam::translateDNA(temp$codon) %in% c("F", "W")), ]
      }
    } else if(is_v_con){
      # V conserved C site landed in J region: allow reference codon OR C
      if(length(j_groups[[i]]) == 3){
        if(sum("N" %in% ref_j[j_groups[[i]]]) == 0){
          if(alakazam::translateDNA(paste0(ref_j[j_groups[[i]]], collapse = "")) != "*"){
            temp <- temp[temp$codon == paste0(ref_j[j_groups[[i]]], collapse = "") |
                           (alakazam::translateDNA(temp$codon) == "C"), ]
          } else{
            temp <- temp[alakazam::translateDNA(temp$codon) == "C", ]
          }
        } else{
          values <- ref_j[j_groups[[i]]]
          values <- paste0(values[-which(values == "N")], collapse = "")
          temp <- temp[(startsWith(temp$codon, values)) |
                         (alakazam::translateDNA(temp$codon) == "C"), ]
        }
      } else{
        values <- paste0(ref_j[j_groups[[i]]], collapse = "")
        temp <- temp[(endsWith(temp$codon, values)) |
                       (alakazam::translateDNA(temp$codon) == "C"), ]
      }
    } else{
      if(length(j_groups[[i]]) == 3){
        if(sum("N" %in% ref_j[j_groups[[i]]]) == 0){
          if(alakazam::translateDNA(paste0(ref_j[j_groups[[i]]], collapse = "")) != "*"){
            temp <- temp[temp$codon == paste0(ref_j[j_groups[[i]]], collapse = ""), ]
          } else{
            warning("A stop codon was detected in the reference for clone ", 
                    cons$clone_id, "which may indicate the improper reference is",
                    " being used")
            writeLines(paste("A stop codon was detected in the reference for clone", 
                             cons$clone_id, "which may indicate the improper",
                             "reference is being used"), 
                       file.path(subDir, "j_annotation_warning.txt"))
          }
        } else{
          if(i != length(j_groups) & pad_length > 0){
            values <- ref_j[j_groups[[i]]]
            values <- paste0(values[-which(values == "N")], collapse = "")
            temp <- temp[startsWith(temp$codon, values), ]
          }
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
  regions <- regions[regions!= "gap"]
  
  germlines_values <- get_starting_junction(tree_df, sub, regions, check_genes_val = TRUE)
  
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
  v <- germlines_values$v
  cdr3 <- germlines_values$cdr3
  j <- germlines_values$j
  temp <- data.frame(v = v, j = j, cdr3 = cdr3)
  return(temp)
}

# gets the tree table from the tree building process
# @param dir      The directory that the tree building method saved to 
# @param subDir   The specific directory for the clone of interest 
# @param clone_ids  The specific clone_id of interest 
# @param repertoire_wide  A logical that indicates if tree building was done on 
#                         a repertoire wide level or not 
getTreeTable <- function(dir, subDir, clone_ids, repertoire_wide){
  if(repertoire_wide){
    tree_df <- suppressWarnings(read.table(file = file.path(dir, "sample", "sample_recon_sample",
                                                            paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")), 
                                           header = F, sep = "\t"))
    file.copy(file.path(dir, "sample", "sample_recon_sample",
                        paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")), 
              file.path(subDir, paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")), 
              overwrite = TRUE)
  } else{
    tree_df <- suppressWarnings(read.table(file = file.path(subDir, "sample", "sample_recon_sample",
                                                            paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")), 
                                           header = F, sep = "\t"))
    file.copy(file.path(subDir, "sample", "sample_recon_sample",
                        paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")), 
              file.path(subDir, paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")), 
              overwrite = TRUE)
  }
  
  colnames(tree_df) = c("site", "codon", "partial_likelihood", "log_likelihood_site",
                        "upper_partial_log_likelihood", "upper_partial_likelihood", "equilibrium")
  tree_df$value <- tree_df$partial_likelihood + log(tree_df$equilibrium)
  return(tree_df)
}

# get the v, j, and junction sequences
# @param tree_df    The df asscoiated with the partial likelihood in tree building 
# @param sub        The clones object for a particualr clone 
# @param regions    The regions of sub -- if sub is run with chain == "HL" this 
#                   vector should only have one chain's worth of regions
# @param check_gene_val   A logical to indicate if this is run for check_gene or not 
get_starting_junction <- function(tree_df, sub, regions, check_genes_val = FALSE){
  if(sub$data[[1]]@phylo_seq == "hlsequence"){
    if(check_genes_val){
      cdr3_index <- (min(which(regions == "cdr3")) - 3):(max(which(regions == "cdr3")) + 3)
    } else{
      heavy_r <- regions[1:nchar(sub$data[[1]]@germline)]
      light_r <- regions[(nchar(sub$data[[1]]@germline) + 1): length(regions)]
      cdr3_index <- (min(which(heavy_r == "cdr3")) - 3):(max(which(heavy_r == "cdr3")) + 3)
    }
  } else{
    cdr3_index <- (min(which(regions == "cdr3")) - 3):(max(which(regions == "cdr3")) + 3)
  }
  
  if(sub$data[[1]]@phylo_seq == "sequence"){
    padding <- length(which(strsplit(substring(sub$data[[1]]@germline, max(cdr3_index) + 1,
                                               length(regions)), "")[[1]] == "N"))
  } else if(sub$data[[1]]@phylo_seq == "lsequence"){
    padding <- length(which(strsplit(substring(sub$data[[1]]@lgermline, max(cdr3_index) + 1,
                                               length(regions)), "")[[1]] == "N"))
  }else{
    if(!check_genes_val){
      padding <- length(which(strsplit(substring(sub$data[[1]]@germline, max(cdr3_index) + 1,
                                                 length(heavy_r)), "")[[1]] == "N"))
      padding_light <- length(which(strsplit(substring(sub$data[[1]]@lgermline, max(which(light_r == "cdr3")) + 1,
                                                       length(light_r)), "")[[1]] == "N"))
    } else{
      if(length(regions) == nchar(sub$data[[1]]@germline)){
        padding <- length(which(strsplit(substring(sub$data[[1]]@germline, max(cdr3_index) + 1,
                                                   length(regions)), "")[[1]] == "N"))
      } else{
        padding <- length(which(strsplit(substring(sub$data[[1]]@lgermline, max(cdr3_index) + 1,
                                                   length(regions)), "")[[1]] == "N"))
      }
    }
  }
  
  tree_seq <- paste0(unlist(lapply(unique(tree_df$site), function(x){
    tmp <- tree_df[tree_df$site == x,]
    sub_indx <- which(tmp$value == max(tmp$value))[1]
    tmp$codon[sub_indx]
  })), collapse = "")
  if(sub$data[[1]]@phylo_seq == "hlsequence" && !check_genes_val){
    nsite_heavy <- nchar(sub$data[[1]]@germline)/3
    tree_df_light <- tree_df[tree_df$site >= nsite_heavy,]
    tree_df_light$site <- tree_df_light$site - min(tree_df_light$site)
    tree_df <- tree_df[tree_df$site < nsite_heavy,]
    tree_seq_light <- substring(tree_seq, nsite_heavy*3 + 1, nchar(tree_seq))
    tree_seq <- substring(tree_seq, 1, nsite_heavy*3)
    cdr3 <- paste0(strsplit(tree_seq, "")[[1]][cdr3_index], collapse = "")
  } else{
    cdr3 <- paste0(strsplit(tree_seq, "")[[1]][cdr3_index], collapse = "")
  }
  test_cdr3 <- strsplit(alakazam::translateDNA(cdr3), "")[[1]]
  if(test_cdr3[1] != "C" || !test_cdr3[length(test_cdr3)] %in% c("F", "W")){
    if(sub$data[[1]]@phylo_seq == "sequence"){
      groupedList <- split(1:length(regions), ceiling(seq_along(regions) / 3))
    } else if(sub$data[[1]]@phylo_seq == "hlsequence"){
      if(!check_genes_val){
        groupedList <- split(1:length(heavy_r), ceiling(seq_along(heavy_r) / 3))
      } else{
        groupedList <- split(1:length(regions), ceiling(seq_along(regions) / 3))
      }
    } else if(sub$data[[1]]@phylo_seq == "lsequence"){
      groupedList <- split(1:length(regions), ceiling(seq_along(regions) / 3))
    }
    if(test_cdr3[1] != "C"){
      codon_site <- which(sapply(groupedList, function(group) min(cdr3_index) %in% group))
      sub_tree_df <- dplyr::filter(tree_df, !!rlang::sym("site") == codon_site - 1)
      sub_tree_df$aa <- alakazam::translateDNA(sub_tree_df$codon)
      sub_tree_df <- dplyr::filter(sub_tree_df, !!rlang::sym("aa") == "C")
      sub_tree_df$value <- sub_tree_df$partial_likelihood + log(sub_tree_df$equilibrium)
      value <- sub_tree_df$codon[sub_tree_df$value == max(sub_tree_df$value)]
      cdr3 <- paste0(value[1], substring(cdr3, 4, nchar(cdr3)))
    }
    if(!test_cdr3[length(test_cdr3)] %in% c("F", "W")){
      codon_site <- which(sapply(groupedList, function(group) max(cdr3_index) %in% group))
      sub_tree_df <- dplyr::filter(tree_df, !!rlang::sym("site") == codon_site - 1)
      sub_tree_df$aa <- alakazam::translateDNA(sub_tree_df$codon)
      sub_tree_df <- dplyr::filter(sub_tree_df, !!rlang::sym("aa") %in% c("W", "F"))
      sub_tree_df$value <- sub_tree_df$partial_likelihood + log(sub_tree_df$equilibrium)
      value <- sub_tree_df$codon[sub_tree_df$value == max(sub_tree_df$value)]
      cdr3 <- paste0(substring(cdr3, 1, nchar(cdr3)-3), value[1])
    }
  }
  if(sub$data[[1]]@phylo_seq == "sequence"){
    v_len <- min(cdr3_index)-1
    v <- substring(tree_seq, 1, v_len)
    j_start <- nchar(paste0(v, cdr3, collapse = "")) + 1
    j <- substring(tree_seq, j_start, nchar(tree_seq))
  } else if(sub$data[[1]]@phylo_seq == "hlsequence"){
    v_len <- min(cdr3_index)-1
    v <- substring(tree_seq, 1, v_len)
    j_start <- nchar(paste0(v, cdr3, collapse = "")) + 1
    j <- substring(tree_seq, j_start, nchar(tree_seq))
  } else if(sub$data[[1]]@phylo_seq == "lsequence"){
    v_len <- min(cdr3_index)-1
    v <- substring(tree_seq, 1, v_len)
    j_start <- nchar(paste0(v, cdr3, collapse = "")) +1
    j <- substring(tree_seq, j_start, nchar(tree_seq))
  }
  if(padding > 0){
    j <- substring(j, 1, nchar(j) - padding)
    j <- paste0(j, paste(rep("N", padding), collapse = ""))
  }
  if(sub$data[[1]]@phylo_seq == "hlsequence" && !check_genes_val){
    v_light <- substring(tree_seq_light, 1, sum(light_r %in% c("cdr1", "cdr2", "fwr1", "fwr2", "fwr3")))
    light_cdr3 <- substring(tree_seq_light, nchar(v_light) + 1, nchar(v_light) + sum(light_r == "cdr3"))
    j_light <- substring(tree_seq_light, nchar(v_light) + nchar(light_cdr3) + 1, nchar(tree_seq_light))
    
    last_v_codon <- substring(v_light, nchar(v_light)-2, nchar(v_light))
    first_j_codon <- substring(j_light, 1, 3)
    light_cdr3 <- paste0(last_v_codon, light_cdr3, first_j_codon)
    v_light <- substring(v_light, 1, nchar(v_light)-3)
    j_light <- substring(j_light, 4, nchar(j_light))
    
    if(padding_light > 0){
      j_light <- substring(j_light, 1, nchar(j_light) - padding_light)
      j_light <- paste0(j_light, paste(rep("N", padding_light), collapse = ""))
    }
    
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
  if(sub$data[[1]]@phylo_seq == "hlsequence" && !check_genes_val){
    return(list(v = v, j = j, cdr3 = cdr3, v_light = v_light, 
                j_light = j_light, light_cdr3 = light_cdr3))
  } else{
    return(list(v = v, j = j, cdr3 = cdr3))
  }
}

# helper for finding germline likelihoods
# @param has_mult A df that is the output of createAllGermlines
# @param sub_data_ref An airr table for a specific clone
# @param v_len    The number of characters in V
# @param j_len    The number of characters in J
# @param tree_df_ref  The df asscoiated with the partial likelihood in tree building 
# @param heavy    A logical to indicate if the process should be run on IGH or nonIGH BCRs
# @param chain    The chain option used in formatClones ("H", "HL", or "L")
# @param locus    The name of the locus column in sub_data
# @param clone    The name of the clone_id column in sub_data
# @param split_light  Split by light chain subgroup?
compute_germlines <- function(has_mult, sub_data_ref, v_len, j_len, tree_df_ref, heavy = TRUE, 
                              chain = "H", locus = "locus", germ = "germline_alignment",
                              clone = "clone_id", split_light = FALSE, ...) {
  do.call(rbind, lapply(1:nrow(has_mult), function(z){
    germline_str <- has_mult$germline_d_mask[z]
    if(heavy){
      sub_data_ref$germline_alignment_d_mask[sub_data_ref[[locus]] == "IGH"] <- germline_str
    } else{
      sub_data_ref$germline_alignment_d_mask[sub_data_ref[[locus]] != "IGH"] <- germline_str
    }
    
    sub_ref <- formatClones(sub_data_ref, chain = chain, clone = clone, minseq = 1, 
                            split_light = split_light, ...)
    
    if(heavy){
      base_germ <- sub_ref$data[[1]]@germline
    } else{
      base_germ <- sub_ref$data[[1]]@lgermline
    }
    
    value_v <- substring(base_germ, 1, v_len)
    value_j <- substring(base_germ, nchar(base_germ)-j_len+1, nchar(base_germ))
    
    v_boundary_codon <- ceiling(nchar(value_v) / 3) - 1
    j_start_codon <- max(tree_df_ref$site) - ceiling(nchar(value_j)/3) + 1
    all_sites <- unique(c(0:v_boundary_codon, j_start_codon:max(tree_df_ref$site)))
    
    sub_df <- dplyr::filter(tree_df_ref, !!rlang::sym("site") %in% all_sites) |>
      dplyr::distinct()
    
    # add Ns to the seqs if needed
    if(nchar(value_v) %% 3 != 0){
      diff <- 3 - nchar(value_v) %% 3
      value_v <- paste0(value_v, paste(rep("N", diff), collapse = ""))
    }
    if(nchar(value_j) %% 3 != 0){
      diff <- 3 - nchar(value_j) %% 3
      value_j <- paste0(paste(rep("N", diff), collapse = ""), value_j)
    }
    vj <- paste0(value_v, value_j)
    gene_list <- strsplit(vj, "")[[1]]
    groupedSeq <- split(gene_list, ceiling(seq_along(gene_list) / 3))
    
    stopifnot(length(groupedSeq) == length(all_sites))
    
    contains_N <- sapply(groupedSeq, function(x) any(x == "N"))
    if(any(contains_N)){
      for(i in which(contains_N) - 1){
        pattern <- groupedSeq[[i + 1]]
        non_n_pos <- which(pattern != "N")
        df_row <- sub_df[sub_df$site == i, ]
        condition <- rep(TRUE, nrow(df_row))
        for(pos in non_n_pos)
          condition <- condition & (substr(df_row$codon, pos, pos) == pattern[pos])
        value_i <- sum(df_row[condition, ]$value)
        new_row <- sub_df[1, ]
        new_row$site <- i
        new_row$codon <- paste0(pattern, collapse = "")
        new_row$value <- value_i
        sub_df <- rbind(sub_df, new_row)
      }
    }
    
    likelihood <- sum(unlist(lapply(1:length(groupedSeq), function(i){
      codon   <- paste0(groupedSeq[[i]], collapse = "")
      sitedf  <- sub_df[sub_df$site == all_sites[i], ]
      sitedf$value[sitedf$codon == codon]
    })))
    
    # just in case 
    value_v <- sub("N+$", "", value_v)
    value_j <- sub("^N+", "", value_j)
    
    data.frame(clone_id  = z, likelihood = likelihood,
               v = value_v, j = value_j, v_call = has_mult$v_call[z], 
               j_call = has_mult$j_call[z], v_start = has_mult$v_start[z], 
               v_end = has_mult$v_end[z], j_start = has_mult$j_start[z],   
               j_end = has_mult$j_end[z], germline = has_mult$germline[z])
  }))
}

# helper pad germline d mask and ungapped to boundary 
# @param vec the string that is to be padded
pad_to_codon <- function(vec){
  sapply(vec, function(x){
    rem <- nchar(x) %% 3
    if(rem != 0) paste0(x, paste(rep("N", 3 - rem), collapse = "")) else x
  })
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
# @param resolve_germ Resolve the V and J genes within the clone?
# @param all_germlines The germlines table needed for resolve_germ
# @param v_call     The name of the v annotation column
# @param j_call     The name of the j annotation column 
# @param locus      The name of the locus column 
# @param quiet      How much noise to print out
# @param clone      The name of the proper clone_id column to use
# @param data_clone The name of the clone id in data that matches clones$clone_id (either the same as clone or clone_subgroup_id)
# @param seq_id     The name of the sequence_id column
# @param chain      H, HL, or L?
# @param check_genes Check if the inferred V/J lengths go into the inferred cdr3 region and adjust accordingly. 
# @param references IMGT references read in using \link{readIMGT}
# @param repertoire_wide Were the trees made using repertoire_wide parameters or do they need to made still?
# @param exec       The exec file path for the tree building method (igphyml)
# @param partition The partition to use when building the tree
# @param v_germ_start V germline start column 
# @param v_germ_end   V germline end column 
# @param j_germ_start J germline start column 
# @param j_germ_end   J germline end column 
# @param germ_align   The germline alignment column 
# @param germ_mask    The gemrline_d_mask column 
# @param search       method for the OLGA bit of the UCA process
# @param cell         The cell id column name 
# @param split_light  A logical to indicate if a clone shoule be split by light chain subgroup or not
processCloneGermline <- function(clone_ids, clones, data, dir, id,
                                 resolve_germ = FALSE, all_germlines = NULL, 
                                 v_call = "v_call", j_call = "j_call", 
                                 locus = "locus", quiet = 0, clone = "clone_id",
                                 data_clone = "clone_id", seq_id = "sequence_id", 
                                 chain = "H", check_genes = FALSE, 
                                 repertoire_wide = FALSE, exec = NULL, 
                                 partition = "single", 
                                 v_germ_start = "v_germline_start", 
                                 v_germ_end = "v_germline_end", 
                                 j_germ_start = "j_germline_start", 
                                 j_germ_end = "j_germline_end", 
                                 germ_align = "germline_alignment",
                                 germ_mask = "germline_alignment_d_mask",
                                 search = "codon", cell = "cell_id", 
                                 split_light = FALSE, ...){
  sub <- dplyr::filter(clones, !!rlang::sym("clone_id") == clone_ids)
  sub_data <- dplyr::filter(data, !!rlang::sym(data_clone) == clone_ids)
  subDir <- file.path(dir, paste0(id, "_",clone_ids))
  if(!dir.exists(subDir)){
    dir.create(subDir, recursive = T)
  }
  
  if(quiet > 0){
    print(paste("constructing trees for", clone_ids))
  }
  
  if(!repertoire_wide){
    sub <- tryCatch({
      getTrees(sub, build = "igphyml", exec = exec, rm_temp = FALSE, dir = subDir,
               asrp = TRUE, nproc = 1, partition = partition, ...)
    }, error = function(e){
      message(paste("getTrees failed for clone", clone_ids, "--retrying,",
                    "Error was:", conditionMessage(e)))
      tryCatch({
        getTrees(sub, build = "igphyml", exec = exec, rm_temp = FALSE, dir = subDir,
                 asrp = TRUE, nproc = 1, partition = partition, ...)
      }, error = function(e2){
        message(paste("getTrees failed twice for clone", clone_ids, "--skipping.", 
                      "Final error was:", conditionMessage(e2)))
        return(NULL)
      })
    })
    if(is.null(sub)){
      return(NULL)
    }
  }
  
  saveRDS(sub, file.path(subDir, "clone.rds"))
  
  if(sub$data[[1]]@phylo_seq == "hlsequence"){
    test_hl <- paste0(strsplit(
      sub$data[[1]]@hlgermline, "")[[1]][(nchar(sub$data[[1]]@germline) + 1):
                                           nchar(sub$data[[1]]@hlgermline)],
      collapse = "")
    if(nchar(test_hl) > nchar(sub$data[[1]]@lgermline)){
      sub$data[[1]]@lgermline <- test_hl
    }
  }
  
  if(sub$data[[1]]@phylo_seq == "hlsequence"){
    heavy_sub <- sub_data[sub_data[[locus]] == "IGH",]
    light_sub <- sub_data[sub_data[[locus]] != "IGH",]
    
    if (nrow(heavy_sub) == 0){
      stop("Clone ", clone_ids, " has no IGH sequences in sub_data. ",
           "Check locus column values.")
    }
    if (nrow(light_sub) == 0){
      stop("Clone ", clone_ids, " has no light-chain sequences in sub_data. ",
           "Check locus column values.")
    }
    
    cons <- findConsensus(heavy_sub, v_call = v_call, j_call = j_call, id = seq_id, ...)
    cons <- sub_data[sub_data[[seq_id]] == cons$cons_id,] 
    cons_light <- findConsensus(light_sub, v_call = v_call, j_call = j_call, id = seq_id, ...)
    cons_light <- sub_data[sub_data[[seq_id]] == cons_light$cons_id,] 
    
  } else if(sub$data[[1]]@phylo_seq == "sequence"){
    cons <- findConsensus(sub_data[sub_data[[locus]] == "IGH",], v_call = v_call,
                          j_call = j_call, id = seq_id, ...)
    cons <- sub_data[sub_data[[seq_id]] == cons$cons_id,]
  } else{
    cons <- findConsensus(sub_data[sub_data[[locus]] != "IGH",], v_call = v_call,
                          j_call = j_call, id = seq_id, ...)
    cons <- sub_data[sub_data[[seq_id]] == cons$cons_id,]
  }
  
  regions <- sub$data[[1]]@region
  
  if(sub$data[[1]]@phylo_seq == "sequence"){
    cdr3_index <- (min(which(regions == "cdr3")) - 3):(max(which(regions == "cdr3")) + 3)
  } else if(sub$data[[1]]@phylo_seq == "hlsequence"){
    heavy_r <- regions[1:nchar(sub$data[[1]]@germline)]
    light_r <- regions[(nchar(sub$data[[1]]@germline) + 1): length(regions)]
    cdr3_index <- (min(which(heavy_r == "cdr3")) - 3):(max(which(heavy_r == "cdr3")) + 3)
  } else if(sub$data[[1]]@phylo_seq == "lsequence"){
    cdr3_index <- (min(which(regions == "cdr3")) - 3):(max(which(regions == "cdr3")) + 3)
  }
  
  tree_df <- getTreeTable(dir, subDir, clone_ids, repertoire_wide)
  
  germline_values <- get_starting_junction(tree_df, sub, regions)
  v <- germline_values$v
  j <- germline_values$j
  cdr3 <- germline_values$cdr3
  
  if(sub$data[[1]]@phylo_seq == "hlsequence"){
    v_light <- germline_values$v_light
    j_light <- germline_values$j_light
    light_cdr3 <- germline_values$light_cdr3
    nsite_heavy <- nchar(sub$data[[1]]@germline)/3
    tree_df_light <- tree_df[tree_df$site >= nsite_heavy,]
    tree_df <- tree_df[tree_df$site < nsite_heavy,]
    tree_df_light$site <- tree_df_light$site - min(tree_df_light$site)
  }
  
  if(resolve_germ){
    if(quiet > 0){
      print(paste("resolving genes for", sub$clone_id))
    }
    
    all_germ_clone <- all_germlines[all_germlines$clone_id == sub$clone_id,]
    heavy_indx <- grepl("^IGH", all_germ_clone$v_call)
    saveRDS(all_germ_clone, file.path(subDir, "all_germlines.rds"))
    
    needs_resolve <- switch(sub$data[[1]]@phylo_seq,
                            "hlsequence" = nrow(all_germ_clone[ heavy_indx, ]) > 1 ||
                              nrow(all_germ_clone[!heavy_indx, ]) > 1,
                            "lsequence" = nrow(all_germ_clone[!heavy_indx, ]) > 1,
                            nrow(all_germ_clone[heavy_indx, ]) > 1)
    
    if(!needs_resolve){
      if(quiet > 0){
        print(paste("only one germline option for", sub$clone_id, "-- skipping resolution"))
      }
      if(sub$data[[1]]@phylo_seq == "hlsequence"){
        has_multiple_light <- all_germ_clone[!heavy_indx,]
        has_multiple_light$ungapped <- pad_to_codon(has_multiple_light$ungapped)
        germlines_light <- compute_germlines(has_mult = has_multiple_light, 
                                             sub_data_ref = sub_data,
                                             v_len = nchar(v_light), 
                                             j_len = nchar(j_light),
                                             tree_df_ref = tree_df_light, 
                                             heavy = FALSE, chain = chain, 
                                             locus = locus, germ = germ_align,
                                             split_light = split_light, ...)
        has_multiple_light$germline_d_mask <- pad_to_codon(has_multiple_light$germline_d_mask)
        saveRDS(germlines_light, file.path(subDir, "most_likely_germlines_light.rds"))
      } 
      
      if(chain == "L"){
        heavy_val = FALSE
      } else{
        heavy_val = TRUE
      }
      
      has_multiple <- switch(sub$data[[1]]@phylo_seq,
                             "hlsequence" = all_germ_clone[heavy_indx, ],
                             "lsequence"  = all_germ_clone[!heavy_indx, ],
                             all_germ_clone[heavy_indx, ])
      
      has_multiple$germline_d_mask <- pad_to_codon(has_multiple$germline_d_mask)
      has_multiple$ungapped <- pad_to_codon(has_multiple$ungapped)
      germlines <- compute_germlines(has_mult = has_multiple, 
                                     sub_data_ref = sub_data,
                                     v_len = nchar(v), j_len = nchar(j), 
                                     tree_df_ref = tree_df, heavy = heavy_val, 
                                     chain = chain, locus = locus, 
                                     germ = germ_align, 
                                     split_light = split_light, ...)
      saveRDS(germlines, file.path(subDir, "most_likely_germlines.rds"))
    } else{
      if(sub$data[[1]]@phylo_seq == "hlsequence"){
        has_multiple_light <- all_germ_clone[!heavy_indx,]
        has_multiple_light$ungapped <- pad_to_codon(has_multiple_light$ungapped)
        has_multiple_light$germline_d_mask <- pad_to_codon(has_multiple_light$germline_d_mask)
        saveRDS(has_multiple_light, file.path(subDir, "all_germlines_light.rds"))
        
        germlines_light <- compute_germlines(has_mult = has_multiple_light, 
                                             sub_data_ref = sub_data,
                                             v_len = nchar(v_light), 
                                             j_len = nchar(j_light),
                                             tree_df_ref = tree_df_light, 
                                             heavy = FALSE, chain = chain, 
                                             locus = locus, germ = germ_align,
                                             split_light = split_light, ...)
        
        index <- which.max(germlines_light$likelihood)
        germlines_light <- germlines_light[index, ]
        v_light <- substring(germlines_light$v, 1, nchar(v_light))
        j_light <- substring(germlines_light$j, nchar(germlines_light$j) - nchar(j_light) + 1,
                             nchar(germlines_light$j))
        saveRDS(germlines_light, file.path(subDir, "most_likely_germlines_light.rds"))
      } else{
        v_light <- NULL
        j_light <- NULL
      }
      
      has_multiple <- switch(sub$data[[1]]@phylo_seq,
                             "hlsequence" = all_germ_clone[heavy_indx, ],
                             "lsequence"  = all_germ_clone[!heavy_indx, ],
                             all_germ_clone[heavy_indx, ])
      
      has_multiple$germline_d_mask <- pad_to_codon(has_multiple$germline_d_mask)
      has_multiple$ungapped <- pad_to_codon(has_multiple$ungapped)
      saveRDS(has_multiple, file.path(subDir, "all_germlines.rds"))
      
      if(chain == "L"){
        heavy_val = FALSE
      } else{
        heavy_val = TRUE
      }
      
      germlines <- compute_germlines(has_mult = has_multiple, 
                                     sub_data_ref = sub_data,
                                     v_len = nchar(v), j_len = nchar(j), 
                                     tree_df_ref = tree_df, heavy = heavy_val, 
                                     chain = chain, locus = locus, 
                                     germ = germ_align, 
                                     split_light = split_light, ...)
      
      index <- which.max(germlines$likelihood)
      germlines <- germlines[index, ]
      v <- germlines$v
      j <- germlines$j
      saveRDS(germlines, file.path(subDir, "most_likely_germlines.rds"))
      
      if(sub$data[[1]]@phylo_seq == "sequence" || sub$data[[1]]@phylo_seq == "hlsequence"){
        cons[[v_call]] <- germlines$v_call
        cons[[j_call]] <- germlines$j_call
        sub_data[[v_call]][sub_data[[locus]] == "IGH"] <- germlines$v_call
        sub_data[[v_germ_start]][sub_data[[locus]] == "IGH"] <- germlines$v_start
        sub_data[[v_germ_end]][sub_data[[locus]] == "IGH" ] <- germlines$v_end
        sub_data[[j_call]][sub_data[[locus]] == "IGH"] <- germlines$j_call
        sub_data[[j_germ_start]][sub_data[[locus]] == "IGH"] <- germlines$j_start
        sub_data[[j_germ_end]][sub_data[[locus]] == "IGH"] <- germlines$j_end
        indx <- which(has_multiple$v_call == germlines$v_call & 
                        has_multiple$j_call == germlines$j_call)
        sub_data[[germ_align]][sub_data[[locus]] == "IGH"] <- has_multiple$germline[indx]
        sub_data[[germ_mask]][sub_data[[locus]] == "IGH"] <- has_multiple$germline_d_mask[indx]
      } else if(sub$data[[1]]@phylo_seq == "lsequence"){
        cons[[v_call]] <- germlines$v_call
        cons[[j_call]] <- germlines$j_call
        sub_data[[v_call]][sub_data[[locus]] != "IGH"] <- germlines$v_call
        sub_data[[v_germ_start]][sub_data[[locus]] != "IGH"] <- germlines$v_start
        sub_data[[v_germ_end]][sub_data[[locus]] != "IGH"] <- germlines$v_end
        sub_data[[j_call]][sub_data[[locus]] != "IGH"] <- germlines$j_call
        sub_data[[j_germ_start]][sub_data[[locus]] != "IGH"] <- germlines$j_start
        sub_data[[j_germ_end]][sub_data[[locus]] != "IGH"] <- germlines$j_end
        indx <- which(has_multiple$v_call == germlines$v_call & 
                        has_multiple$j_call == germlines$j_call)
        sub_data[[germ_align]][sub_data[[locus]] != "IGH"] <- has_multiple$germline[indx]
        sub_data[[germ_mask]][sub_data[[locus]] != "IGH"] <- has_multiple$germline_d_mask[indx]
      }
      
      if(sub$data[[1]]@phylo_seq == "hlsequence"){
        cons_light[[v_call]] <- germlines_light$v_call
        cons_light[[j_call]] <- germlines_light$j_call
        sub_data[[v_call]][sub_data[[locus]] != "IGH"] <- germlines_light$v_call
        sub_data[[v_germ_start]][sub_data[[locus]] != "IGH"] <- germlines_light$v_start
        sub_data[[v_germ_end]][sub_data[[locus]] != "IGH"] <- germlines_light$v_end
        sub_data[[j_call]][sub_data[[locus]] != "IGH"] <- germlines_light$j_call
        sub_data[[j_germ_start]][sub_data[[locus]] != "IGH"] <- germlines_light$j_start
        sub_data[[j_germ_end]][sub_data[[locus]] != "IGH"] <- germlines_light$j_end
        indx <- which(has_multiple_light$v_call == germlines_light$v_call & 
                        has_multiple_light$j_call == germlines_light$j_call)
        sub_data[[germ_align]][sub_data[[locus]] != "IGH"] <- has_multiple_light$germline[indx]
        sub_data[[germ_mask]][sub_data[[locus]] != "IGH"] <- has_multiple_light$germline_d_mask[indx]
      }
    }

    if(quiet > 0){
      print("reconstructing trees")
    }
    
    saveRDS(sub, file.path(subDir, "original_clone.rds"))
    sub_ids <- sub$data[[1]]@data$sequence_id
    sub_data_clone <- sub_data[sub_data[[seq_id]] %in% sub_ids,]
    if(cell %in% colnames(sub_data)){
      sub_data_clone <- sub_data[sub_data[[cell]] %in% sub_data_clone[[cell]],]
    }
    
    sub <- formatClones(sub_data_clone, chain = chain, clone = clone, 
                        dup_singles = TRUE, minseq = 1, 
                        split_light = split_light, ...)
    
    # rename the old folder 
    file.rename(file.path(subDir, "sample"), file.path(subDir, "masked_sample"))
    
    sub <- tryCatch({
      getTrees(sub, build = "igphyml", exec = exec, rm_temp = FALSE, dir = subDir,
               asrp = TRUE, nproc = 1, partition = partition, ...)
    }, error = function(e){
      message(paste("getTrees failed for clone", clone_ids, "--retrying,",
                    "Error was:", conditionMessage(e)))
      tryCatch({
        getTrees(sub, build = "igphyml", exec = exec, rm_temp = FALSE, dir = subDir,
                 asrp = TRUE, nproc = 1, partition = partition, ...)
      }, error = function(e2){
        message(paste("getTrees failed twice for clone", clone_ids, "--skipping.", 
                      "Final error was:", conditionMessage(e2)))
        return(NULL)
      })
    })
      
    tree_df <- getTreeTable(dir, subDir, clone_ids, repertoire_wide)
    regions <- sub$data[[1]]@region
    
    if(sub$data[[1]]@phylo_seq == "hlsequence"){
      heavy_r <- regions[1:nchar(sub$data[[1]]@germline)]
      light_r <- regions[(nchar(sub$data[[1]]@germline) + 1): length(regions)]
      cdr3_index <- (min(which(heavy_r == "cdr3")) - 3):(max(which(heavy_r == "cdr3")) + 3)
    }
    
    germline_values <- get_starting_junction(tree_df, sub, regions)
    v <- germline_values$v
    j <- germline_values$j
    cdr3 <- germline_values$cdr3
    
    if(sub$data[[1]]@phylo_seq == "hlsequence"){
      v_light <- germline_values$v_light
      j_light <- germline_values$j_light
      light_cdr3 <- germline_values$light_cdr3
      nsite_heavy <- nchar(sub$data[[1]]@germline)/3
      tree_df_light <- tree_df[tree_df$site >= nsite_heavy,]
      tree_df <- tree_df[tree_df$site < nsite_heavy,]
      tree_df_light$site <- tree_df_light$site - min(tree_df_light$site)
    } 
  }
  
  saveRDS(sub, file.path(subDir, "clone.rds"))
  
  if(check_genes){
    if(quiet > 0){
      print(paste("restricting tree options to germline references for", sub$clone_id))
    }
    if(sub$data[[1]]@phylo_seq == "hlsequence"){
      heavy_vals <- checkGenesUCA(sub = sub, cons = cons,  v = v, cdr3 = cdr3, 
                                  j = j, tree_df = tree_df, subDir = subDir,
                                  clone_ids = clone_ids, chain = "H", 
                                  regions = heavy_r, clone = clone, seq_id = seq_id)
      light_vals <- checkGenesUCA(sub = sub, cons = cons_light, v = v_light, cdr3 = light_cdr3,
                                  j = j_light, tree_df = tree_df_light, subDir = subDir, 
                                  clone_ids = clone_ids, chain = "L", 
                                  regions = light_r, clone = clone, seq_id = seq_id)
      v <- heavy_vals$v
      cdr3 <- heavy_vals$cdr3
      j <- heavy_vals$j
      v_light <- light_vals$v
      light_cdr3 <- light_vals$cdr3
      j_light <- light_vals$j
    } else if(sub$data[[1]]@phylo_seq == "sequence"){
      heavy_vals <- checkGenesUCA(sub = sub, cons = cons, v = v, cdr3 = cdr3, 
                                  j = j, tree_df = tree_df, subDir = subDir, 
                                  clone_ids = clone_ids, chain = "H",
                                  regions = regions, clone = clone, seq_id = seq_id)
      v <- heavy_vals$v
      cdr3 <- heavy_vals$cdr3
      j <- heavy_vals$j
    } else if(sub$data[[1]]@phylo_seq == "lsequence"){
      light_vals <- checkGenesUCA(sub = sub, cons = cons, v = v, cdr3 = cdr3, 
                                  j = j, tree_df = tree_df, subDir = subDir, 
                                  clone_ids = clone_ids, chain = "L",
                                  regions = regions, clone = clone, seq_id = seq_id)
      v <- light_vals$v
      cdr3 <- light_vals$cdr3
      j <- light_vals$j
    }
  }
  
  if(quiet > 0){
    print(paste("sucessfully obtained most likely junction for", clone_ids))
  }
  
  # put it all together 
  v_cdr3 <- paste0(v, paste0(cdr3, collapse = ""), collapse = "")
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
callOlga <- function(clones, dir, uca_script, python, max_iters, nproc, id, model_folder,
                     model_folder_igk, model_folder_igl, quiet, search, ...){
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
  
  args_keys <- args[seq(1, length(args), 2)]
  args_values <- args[seq(2, length(args), 2)]
  args_list <- stats::setNames(as.list(args_values), gsub("^--", "", args_keys))
  
  json_str <- paste0("{\n", paste(sprintf(
    '  "%s": "%s"', names(args_list), gsub('"', '\\"', args_list)), collapse = ",\n"),
    "\n}")
  
  writeLines(json_str, file.path(path.expand(dir), "olga_args.json"))
  
  cmd <- paste(
    shQuote(python), 
    shQuote(path.expand(uca_script)), 
    "--args_json", shQuote(file.path(path.expand(dir), "olga_args.json")))
  
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
# @param data       The AIRR-table associated with the clones object
# @param references The output of readIMGT() -- germline references
# @param dir        The directory where data should be saved to 
# @param id         The run id
# @param nproc      The number of cores to use 
# @param fill_gaps  A logical to fill in the UCA sequences with partial V/J genes. 
# @param clone_id   The clone id varaible 
# @param locus      The locus column 
# @param seq_id     The sequence id 
# @param resolve_germ     A logical for it resolve germ was run 
# @param v_start    Where the v germline starts
# @param v_end      Where the v germline end
# @param j_start    Where the j germline starts
# @param j_end      Where the j germline end
# @param v_call     The v call 
# @param j_call     The j call 
updateClone <- function(clones, data, references, dir, id, nproc = 1,
                        fill_gaps = TRUE, clone_id = "clone_id", locus = "locus",
                        seq_id = "sequence_id", resolve_germ = FALSE,
                        v_start = "v_germline_start", v_end = "v_germline_end",
                        j_start = "j_germline_start", j_end = "j_germline_end",
                        v_call = "v_call", j_call = "j_call"){
  updated_clones <- do.call(rbind, parallel::mclapply(1:nrow(clones), function(x){
    clone <- clones[x,]
    # get locus value from clone -- use as name for uca/uca_aa
    locus_val <- strsplit(clone$locus, ",")[[1]]
    if(clone$data[[1]]@phylo_seq == "lsequence"){
      uca <- read.table(file.path(dir, paste0(id, "_", clone$clone_id),
                                  "UCA_light.txt"), sep = "\t")[[1]]
      names(uca) <- locus_val
    } else{
      uca <- read.table(file.path(dir, paste0(id, "_", clone$clone_id), "UCA.txt"), sep = "\t")[[1]]
      names(uca) <- locus_val[1]
      if(clone$data[[1]]@phylo_seq == "hlsequence"){
        uca_light <- read.table(file.path(dir, paste0(id, "_", clone$clone_id),
                                          "UCA_light.txt"), sep = "\t")[[1]]
        names(uca_light) <- locus_val[2]
        uca <- append(uca, uca_light)
      }
    }
    
    uca_aa <- alakazam::translateDNA(uca)
    
    germline_node <- ape::getMRCA(clone$trees[[1]], clone$trees[[1]]$tip.label)
    if(clone$data[[1]]@phylo_seq != "hlsequence"){
      clone$trees[[1]]$nodes[[germline_node]]$sequence <- uca
    } else{
      clone$trees[[1]]$nodes[[germline_node]]$sequence <- paste0(uca[1], uca[2])
    }
    
    if(clone$data[[1]]@phylo_seq == "hlsequence"){
      UCA_heavy <- getNodeSeq(clone, node = germline_node, tree = clone$trees[[1]])[1]
      UCA_light <- getNodeSeq(clone, node = germline_node, tree = clone$trees[[1]])[2]
      UCA_h_aa <- alakazam::translateDNA(UCA_heavy)
      UCA_l_aa <- alakazam::translateDNA(UCA_light)
      uca_gapped <- c(UCA_heavy, UCA_light)
      uca_gapped_aa <- alakazam::translateDNA(uca_gapped)
    } else if(clone$data[[1]]@phylo_seq %in% c("sequence", "lsequence")){
      uca_gapped <- getNodeSeq(clone, node = germline_node, tree = clone$trees[[1]])[1]
      uca_gapped_aa <- alakazam::translateDNA(uca_gapped)
    } 
    
    if(fill_gaps){
      full_ucas <- c()
      for(i in 1:length(uca_gapped)){
        locus_val <- names(uca_gapped[i])
        c_uca <- strsplit(uca_gapped[i], "")[[1]]
        germ_values <- findConsensus(data[data[[clone_id]] == clone$clone_id &
                                            data[[locus]] == locus_val, ])
        cell_df <- data[data[[seq_id]] == germ_values$cons_id, ]
        if(resolve_germ){
          if(locus_val == "IGH" | clone$data[[1]]@phylo_seq == "lsequence"){
            ml_germ <- readRDS(file.path(dir, paste0(id, "_", clone$clone_id), "most_likely_germlines.rds"))
          } else{
            ml_germ <- readRDS(file.path(dir, paste0(id, "_", clone$clone_id), "most_likely_germlines_light.rds"))
          }
          v_indx <- which(names(references[[locus_val]]$V) == ml_germ$v_call)
          j_indx <- which(names(references[[locus_val]]$J) == ml_germ$j_call)
          v_ref <- references[[locus_val]]$V[v_indx]
          j_ref <- references[[locus_val]]$J[j_indx]
          v_ref <- substring(v_ref, ml_germ$v_start, ml_germ$v_end)
          j_ref <- substring(j_ref, ml_germ$j_start, ml_germ$j_end)
          v_gaps <- which(strsplit(v_ref, "")[[1]] == ".")
          j_gaps <- which(strsplit(j_ref, "")[[1]] == ".")
          gaps <- append(v_gaps, j_gaps)
          clone_germ <- ml_germ$germline
        } else{
          v_indx <- which(names(references[[locus_val]]$V) == strsplit(cell_df[[v_call]], ",")[[1]][1])
          j_indx <- which(names(references[[locus_val]]$J) == strsplit(cell_df[[j_call]], ",")[[1]][1])
          v_ref <- references[[locus_val]]$V[v_indx]
          j_ref <- references[[locus_val]]$J[j_indx]
          v_ref <- substring(v_ref, cell_df[[v_start]], cell_df[[v_end]])
          j_ref <- substring(j_ref, cell_df[[j_start]], cell_df[[j_end]])
          v_gaps <- which(strsplit(v_ref, "")[[1]] == ".")
          j_gaps <- which(strsplit(j_ref, "")[[1]] == ".")
          gaps <- append(v_gaps, j_gaps)
          clone_germ <- cell_df$germline_alignment
        }
        clone_gaps <- dplyr::setdiff(1:max(clone$data[[1]]@numbers), 
                                     clone$data[[1]]@numbers)
        if(!identical(clone_gaps, gaps)){
          # get the values that are not found in the current UCA
          missing_values <- dplyr::setdiff(clone_gaps, gaps)
          germ_gapped <- strsplit(clone_germ, "")[[1]]
          c_uca[missing_values] <- germ_gapped[missing_values]
          c_uca <- paste0(c_uca, collapse = "")
        } else{
          c_uca <- paste0(c_uca, collapse = "")
        }
        names(c_uca) <- locus_val
        full_ucas <- append(full_ucas, c_uca)
      }
      uca_gapped <- full_ucas
      uca_gapped_aa <- alakazam::translateDNA(uca_gapped)
    } 
    clone$UCA <- list(list(
      ungapped = uca,
      ungapped_aa = uca_aa, 
      gapped = uca_gapped, 
      gapped_aa = uca_gapped_aa
    ))
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
#' @param genotyped     A logical that indicates if the data has been genotyped. 
#'                      Will restrict created germlines to be the genotypes found in entire repertoire. 
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
                               verbose=0, genotyped = FALSE, ...){
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
  if(genotyped){
    all_v <- unique(unlist(lapply(data[[v_call]],function(x)
      alakazam::getAllele(x, strip_d=FALSE, first = TRUE))))
    all_j <- unique(unlist(lapply(data[[j_call]],function(x)
      alakazam::getAllele(x, strip_d=FALSE, first = TRUE))))
    genotype_list <- c(all_v, all_j)
  } else{
    genotype_list <- NULL
  }
  unique_clones <- unique(data[,unique(c(clone,fields)),drop=F])
  germlines <- do.call(rbind, parallel::mclapply(1:nrow(unique_clones), function(x){
    sub <- dplyr::right_join(data, unique_clones[x,,drop=F], by=c(clone,fields))
    if(verbose > 0){
      print(x)
    }
    sub_germlines <- do.call(rbind, lapply(unique(sub[[locus]]), function(l){
      buildAllClonalGermlines(receptors = sub[sub[[locus]] == l,],
                              references = references, genotype_list = genotype_list, chain = l, 
                              seq = seq, v_call = v_call, d_call = d_call, j_call = j_call, 
                              amino_acid = amino_acid, id = id, clone = clone, v_germ_start = v_germ_start,
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
#' @param genotype_list    A list of different V and J calls based on genotyping results 
#' @param ...              Additional arguments passed to \link{buildGermline}
#' @return A data frame with all possible reconstructed germlines
#' @seealso \link{createAllGermlines} \link{buildGermline}, \link{stitchVDJ}
#' @export
#' 
buildAllClonalGermlines <- function(receptors, references, genotype_list = NULL,
                                    chain="IGH", use_regions=FALSE, vonly=FALSE,
                                    seq="sequence_alignment", id="sequence_id", clone="clone_id",
                                    v_call="v_call", j_call="j_call",  v_germ_start="v_germline_start",
                                    v_germ_end="v_germline_end", v_germ_length="v_germline_length", 
                                    d_germ_start="d_germline_start", d_germ_end="d_germline_end", 
                                    d_germ_length="d_germline_length", j_germ_start="j_germline_start", 
                                    j_germ_end="j_germline_end", j_germ_length="j_germline_length", 
                                    np1_length="np1_length", np2_length="np2_length",
                                    j_germ_aa_length="j_germline_aa_length", 
                                    amino_acid=FALSE, ...){
  
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
  
  # remove any padding on the sequence alignment (in case it was already padded)
  receptors[[seq]] <- unlist(lapply(1:nrow(receptors), function(x){
    if(is.na(receptors[[np1_length]][x])){
      receptors[[np1_length]][x] <- 0
    }
    if(is.na(receptors[[d_germ_length]][x])){
      receptors[[d_germ_length]][x] <- 0
    }
    if(is.na(receptors[[np2_length]][x])){
      receptors[[np2_length]][x] <- 0
    }
    germline_len <- sum(receptors[[v_germ_length]][x], receptors[[np1_length]][x],
                        receptors[[d_germ_length]][x], receptors[[np2_length]][x],
                        receptors[[j_germ_length]][x])
    seq_len <- nchar(receptors[[seq]][x])
    if(germline_len != seq_len){
      diff <- abs(germline_len - seq_len)
      value <- substring(receptors[[seq]][x], 1, nchar(receptors[[seq]][x])-diff)
    } else{
      value <- receptors[[seq]][x]
    }
    return(value)
  }))
  
  # has to be first due to the igblast coordinates 
  v_dict <- unlist(lapply(receptors[[v_call]],function(x)
    alakazam::getAllele(x, strip_d=FALSE, first = TRUE))) 
  j_dict <- unlist(lapply(receptors[[j_call]],function(x)
    alakazam::getAllele(x, strip_d=FALSE, first = FALSE)))
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
  
  # info for all the germlines 
  if(is.null(genotype_list)){
    # v_split <- lapply(strsplit(receptors[[v_call]], ",", fixed = TRUE), function(x) sub("\\*.*$", "", x))
    # j_split <- lapply(strsplit(receptors[[j_call]], ",", fixed = TRUE), function(x) sub("\\*.*$", "", x))
    # v_families <- unique(trimws(unlist(v_split, use.names = FALSE)))
    # j_families <- unique(trimws(unlist(j_split, use.names = FALSE)))
    # 
    # # get all the v/j genes within the observed families 
    # v_names <- names(references[[chain]]$V)
    # j_names <- names(references[[chain]]$J)
    # 
    # v_pat <- paste0("^(", paste(v_families, collapse = "|"), ")\\*")
    # j_pat <- paste0("^(", paste(j_families, collapse = "|"), ")\\*")
    # 
    # v_all <- v_names[grepl(v_pat, v_names)]
    # j_all <- j_names[grepl(j_pat, j_names)]
    v_all <- unique(unlist(strsplit(receptors[[v_call]], ",", fixed = TRUE)))
    j_all <- unique(unlist(strsplit(receptors[[j_call]], ",", fixed = TRUE)))
  } else{
    v_indx <- startsWith(genotype_list, paste0(chain, "V"))
    v_all <- genotype_list[v_indx]
    j_all <- genotype_list[!v_indx]
  }
  
  
  # All combinations
  combinations <- expand.grid(v_all, j_all, stringsAsFactors = FALSE)
  colnames(combinations) <- c(v_call, j_call)
  
  all_germlines <- c()
  for(x in 1:nrow(combinations)){
    v_cons <- as.character(combinations[[v_call]][x])
    j_cons <- as.character(combinations[[j_call]][x])
    
    cons <- cons_normal
    cons[[v_call]] <- v_cons
    cons[[j_call]] <- j_cons
    
    
    # make sure the cons has D and NP values 
    if(is.na(cons[[d_germ_length]])) cons[[d_germ_length]] <- 0
    if(is.na(cons[[d_germ_start]])) cons[[d_germ_start]] <- 0
    if(is.na(cons[[d_germ_end]])) cons[[d_germ_end]] <- 0
    if(is.na(cons[[np1_length]])) cons[[np1_length]] <- 0
    if(is.na(cons[[np2_length]])) cons[[np2_length]] <- 0
    
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
    
    if(length(sub_db) == 0){
      stop(paste("Reference database for",chain,"is empty"))
    }
    
    germlines <- tryCatch(buildGermline(cons, references=sub_db, seq=seq, 
                                        v_call=v_call, j_call=j_call, 
                                        v_germ_length = v_germ_length,
                                        np1_length = np1_length,
                                        d_germ_length = d_germ_length,
                                        np2_length = np2_length,
                                        j_germ_length=j_germ_length,
                                        amino_acid=amino_acid),error=function(e)e)
    
    if("error" %in% class(germlines) | "N" %in% strsplit(germlines$vonly, "")[[1]]){
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
                         v_start = cons[[v_germ_start]],
                         v_end = cons[[v_germ_end]], 
                         j_start = cons[[j_germ_start]], 
                         j_end = cons[[j_germ_end]],
                         germline = NA,
                         germline_d_mask = NA, 
                         regions = NA, 
                         positions = NA,
                         npositions = NA)
      all_germlines <- rbind(all_germlines, temp)
      next
    }
    
    positions <- as.numeric(gregexpr("\\.", germlines$full)[[1]])
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
                       v_start = cons[[v_germ_start]],
                       v_end = cons[[v_germ_end]], 
                       j_start = cons[[j_germ_start]], 
                       j_end = cons[[j_germ_end]],
                       germline = germlines$full,
                       germline_d_mask = germlines$dmask, 
                       regions = germlines$regions, 
                       positions = paste0(positions, collapse = ","),
                       npositions = length(positions))
    all_germlines <- rbind(all_germlines, temp)
  }
  
  if(sum(is.na(all_germlines$germline)) > 0){
    all_germlines <- all_germlines[-which(is.na(all_germlines$germline)),]
  }
  
  all_germlines$ungapped <- unlist(lapply(1:nrow(all_germlines), function(x){
    germ <- gsub("\\.", "", all_germlines$germline_d_mask[x])
    return(germ)
  }))
  all_germlines$nchar <- nchar(all_germlines$ungapped)
  indx <- intersect(which(strsplit(cons_normal[[v_call]], ",")[[1]][1] == all_germlines$v_call),
                    which(strsplit(cons_normal[[j_call]], ",")[[1]][1] == all_germlines$j_call))
  base_length <- all_germlines$nchar[indx]
  if(any(abs(all_germlines$nchar - base_length) > 0)){
    all_germlines <- all_germlines[abs(all_germlines$nchar - base_length) == 0,]
  }
  return(all_germlines)
}


# \link{maskAmbigousReferenceSites} Determines and builds all possible germlines for a clone
# @param clones AIRR-table that is the output of \link{formatClones}.
# @param data   The df used to make clones
# @param all_germlines A data frame with all possible reconstructed germlines.
# @param dir  Dir to save to
# @param id   Run id
# @param clone Column name for the clone ID.
# @param locus Column name for the locus in data
# @param split_light A logical to indicate split clones by light chain subgroups or not
# @param nproc Number of cores to use for parallel processing. Default is 1.
# @return A data frame with all possible reconstructed germlines.
maskAmbiguousReferenceSites <- function(clones, data, all_germlines,  
                                       clone = "clone_id", locus = "locus", 
                                       split_light = FALSE, nproc = 1, ...){
  ambig_table <- data.frame(value = c("R", "K", "S", "Y", "M", "W", "B", "H", 
                                      "N", "D", "V"), 
                            combo_l = c(2, 2, 2, 2, 2, 2, 3, 3, 4, 3, 3), 
                            combo = c("A,G", "G,T", "C,G", "C,T", "A,C", 
                                      "A,T", "C,G,T", "A,C,T", "A,C,G,T", 
                                      "A,G,T", "A,C,G"))
  updated_clones <- do.call(rbind, parallel::mclapply(clones$clone_id, function(x){
    sub <- clones[which(clones$clone_id == x),]
    sub_germs <- all_germlines[all_germlines$clone_id == x,]
    sub_data <- data[data[[clone]] == x,]
    
    if(sub$data[[1]]@phylo_seq == "sequence"){
      chain <- "H"
    } else if(sub$data[[1]]@phylo_seq == "hlsequence"){
      chain <- "HL"
    } else if(sub$data[[1]]@phylo_seq == "lsequence"){
      chain <- "L"
    }
    
    heavy_n <- nrow(sub_germs[grepl("^IGH", sub_germs$v_call), ])
    light_n  <- nrow(sub_germs[!grepl("^IGH", sub_germs$v_call), ])
    
    needs_mask <- switch(chain, 
                         "H" = heavy_n > 1,
                         "L" = light_n  > 1,
                         "HL" = heavy_n  > 1 || light_n > 1)
    
    if(!needs_mask) return(sub)
    
    if(chain == "H"){
      heavy_indx <- grepl("^IGH", sub_germs$v_call)
      sub_germs <- sub_germs[heavy_indx,]
      comp_df <- do.call(cbind, lapply(1:nrow(sub_germs), function(z){
        temp <- data.frame(strsplit(sub_germs$germline_d_mask[z], "")[[1]])
      }))
      colnames(comp_df) <- sub_germs$clone_id_unique
    } else if(chain == "L"){
      light_indx <- !grepl("^IGH", sub_germs$v_call)
      sub_germs <- sub_germs[light_indx,]
      comp_df <- do.call(cbind, lapply(1:nrow(sub_germs), function(z){
        temp <- data.frame(strsplit(sub_germs$germline_d_mask[z], "")[[1]])
      }))
      colnames(comp_df) <- sub_germs$clone_id_unique
    } else if(chain == "HL"){
      heavy_indx <- grepl("^IGH", sub_germs$v_call)
      light_indx <- !grepl("^IGH", sub_germs$v_call)
      sub_germs_h <- sub_germs[heavy_indx,]
      sub_germs_l <- sub_germs[light_indx,]
      comp_df_h <- do.call(cbind, lapply(1:nrow(sub_germs_h), function(z){
        temp <- data.frame(strsplit(sub_germs_h$germline_d_mask[z], "")[[1]])
      }))
      colnames(comp_df_h) <- sub_germs_h$clone_id_unique
      
      comp_df_l <- do.call(cbind, lapply(1:nrow(sub_germs_l), function(z){
        temp <- data.frame(strsplit(sub_germs_l$germline_d_mask[z], "")[[1]])
      }))
      colnames(comp_df_l) <- sub_germs_l$clone_id_unique
    }
    
    if(chain %in% c("H", "L")){
      comp_df$diff <- apply(comp_df, 1, function(row) length(unique(row)) > 1)
      masked_sites <- which(comp_df$diff)
      new_clone_germ <- rep(".", nrow(comp_df))
      for(i in 1:nrow(comp_df)){
        row_vals <- comp_df[i,]
        row_vals <- row_vals[,-which(colnames(row_vals) == "diff")]
        unique_vals <- unique(unlist(row_vals))
        if(i %in% masked_sites){
          combo_string <- paste(sort(unique_vals), collapse = ",")
          result <- ambig_table$value[ambig_table$combo == combo_string]
          if(length(result) == 0){
            result <- "N"
          }
        } else{
          result <- unique_vals
        }
        new_clone_germ[i] <- result
      }
      new_clone_germ <- paste0(new_clone_germ, collapse = "")
      if(chain == "H"){
        sub_data$germline_alignment_d_mask[sub_data[[locus]] == "IGH"] <- new_clone_germ
      } else{
        sub_data$germline_alignment_d_mask[sub_data[[locus]] != "IGH"] <- new_clone_germ
      }
    } else{
      #  heavy germline
      comp_df_h$diff <- apply(comp_df_h, 1, function(row) length(unique(row)) > 1)
      masked_sites <- which(comp_df_h$diff)
      new_clone_germ_h <- rep(".", nrow(comp_df_h))
      for(i in 1:nrow(comp_df_h)){
        row_vals <- comp_df_h[i,]
        row_vals <- row_vals[,-which(colnames(row_vals) == "diff")]
        unique_vals <- unique(unlist(row_vals))
        if(i %in% masked_sites){
          combo_string <- paste(sort(unique_vals), collapse = ",")
          result <- ambig_table$value[ambig_table$combo == combo_string]
          if(length(result) == 0){
            result <- "N"
          }
        } else{
          result <- unique_vals
        }
        new_clone_germ_h[i] <- result
      }
      new_clone_germ_h <- paste0(new_clone_germ_h, collapse = "")
      
      #light germline 
      comp_df_l$diff <- apply(comp_df_l, 1, function(row) length(unique(row)) > 1)
      masked_sites <- which(comp_df_l$diff)
      new_clone_germ_l <- rep(".", nrow(comp_df_l))
      for(i in 1:nrow(comp_df_l)){
        row_vals <- comp_df_l[i,]
        row_vals <- row_vals[,-which(colnames(row_vals) == "diff")]
        unique_vals <- unique(unlist(row_vals))
        if(i %in% masked_sites){
          combo_string <- paste(sort(unique_vals), collapse = ",")
          result <- ambig_table$value[ambig_table$combo == combo_string]
          if(length(result) == 0){
            result <- "N"
          }
        } else{
          result <- unique_vals
        }
        new_clone_germ_l[i] <- result
      }
      new_clone_germ_l <- paste0(new_clone_germ_l, collapse = "")
      sub_data$germline_alignment_d_mask[sub_data[[locus]] == "IGH"] <- new_clone_germ_h
      sub_data$germline_alignment_d_mask[sub_data[[locus]] != "IGH"] <- new_clone_germ_l
    }

    sub <- formatClones(sub_data, chain = chain, germ = "germline_alignment_d_mask",
                        clone = clone, nproc = 1, split_light = split_light,
                        minseq = 1, ...)
    
    return(sub)
  }, mc.cores = nproc))
  return(updated_clones)
}

#' \link{getTreesAndUCAs} Construct trees and infer the UCA
#' @param clones        AIRR-table containing sequences \link{formatClones}
#' @param data          The AIRR-table that was used to make the clones object.
#' @param dir           The file path of the directory of where data is saved. NULL is default.
#' @param exec          File path to the tree building executable
#' @param resolve_germ  Resolve the V and J gene annotations within each clone?
#' @param repertoire_wide Build trees using parameters inferred from the entire dataset?
#' @param partition     The partition model to use with IgPhyML. "single" is the default.
#' @param model_folder  The file path to the OLGA default model files for heavy chains
#' @param model_folder_igk  The file path to the OLGA default model files for IGK
#' @param model_folder_igl  The file path to the OLGA default model files for IGL
#' @param python        Specify the python call for your system. This is the call
#'                      on command line that issues the python you want to use. 
#'                      "python3" by default. 
#' @param id            The run ID, sample by default
#' @param max_iters     The maximum number of iterations to run before ending. 
#'                      100 by default
#' @param nproc         The number of cores to use 
#' @param rm_temp       Remove the generated files?
#' @param chain         Set to HL to use both heavy and light chain sequences
#' @param quiet         Amount of noise to print out
#' @param references    Reference genes. See \link{readIMGT}
#' @param clone         The name of the clone id column used in \link{formatClones}. 
#'                      If split_light was used in \link{formatClones}, use "clone_subgroup_id".
#' @param cell          The name of the cell id in the AIRR table used to generate \link{formatClones}
#' @param subsample_size The amount that the clone should be sampled down to. Default is NA. Use NA if you do not wish to subsample -- in testing
#' @param search        Search codon or nt space
#' @param check_genes   Check if the inferred V/J lengths go into the inferred cdr3 region and adjust accordingly. 
#' @param fill_gaps     A logical that will fill in the V and J UCAs of clones that have partial V/J sequence alignments
#' @param genotyped     A logical that indicates that AIRR-table has been genotyped and will only use the found annotations to create various germlines
#' @param split_light   A logical that indicates if different light chain groups should be used to further split a clone (recommended for paired data)
#' @param ...           Additional arguments passed to various other functions like \link{getTrees} and \link{buildGermline}
#' @return An \code{airrClone} object with trees and the inferred UCA
#' @details Return object adds/edits following columns:
#' \itemize{
#'   \item  \code{trees}:  The phylogenies associated with each clone
#'   \item  \code{UCA}:    The inferred UCA
#' }
#' @seealso \link{getTrees} 
#' @export
getTreesAndUCAs <- function(clones, data, dir = NULL, exec = NULL, 
                            resolve_germ = FALSE, repertoire_wide = FALSE, 
                            partition = "single", model_folder, 
                            model_folder_igk = NULL, model_folder_igl = NULL, 
                            python = "python3", id = "sample", 
                            max_iters = 100, nproc = 1, 
                            rm_temp = TRUE, quiet = 0, chain = "H",
                            references = NULL, clone = "clone_id", 
                            cell = "cell_id", subsample_size = NA, 
                            search = "codon", check_genes = TRUE, fill_gaps = TRUE, 
                            genotyped = FALSE, split_light = FALSE, ...){
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
  
  if(resolve_germ & is.null(data) | resolve_germ & is.null(references)){
    stop('resolve_germ requires the data object and references', 
         "References need to be read in using dowser::readIMGT()")
  }
  
  if(chain %in% c("L", "HL") & is.null(model_folder_igk) | chain %in% c("L", "HL") & is.null(model_folder_igl)){
    stop("Light chain model folders are required to run paired or light chain UCA reconstruction")
  }
  
  if(chain %in% c("L", "HL") & !"clone_subgroup_id" %in% colnames(data)){
    stop('To run paired chain analysis you need to first run dowser::resolveLightChains')
  }
  
  if(genotyped){
    warning('Novel alleles are not yet supported')
  }
  
  if(chain == "HL" && !split_light){
    warning(
      paste("Paired trees have been requested but not to split by ligth chain subgroup.", 
            "This will result in only the major subgroup to be constructed."))
  }
  
  if(!cell %in% colnames(data) && chain == "HL" && resolve_germ){
    stop(cell, " is not found in data. Please provide a cell id")
  }
  
  if(split_light){
    # check for clone_subgroup_id
    if(!"clone_subgroup_id" %in% colnames(data)){
      stop('Please run dowser::resolveLightChains and restart')
    }
    data_clone_val <- 'clone_subgroup_id'
  } else{
    data_clone_val <- clone
  }
  
  if(resolve_germ){
    all_germlines <- suppressWarnings(
      createAllGermlines(data = data, references = references, nproc = nproc,
                         clone = data_clone_val, trim_lengths = TRUE, verbose = quiet, 
                         genotyped = genotyped, ...))
    saveRDS(all_germlines, file.path(dir, "all_germlines.rds"))
    clones <- maskAmbiguousReferenceSites(clones = clones, 
                                         all_germlines = all_germlines, data = data,
                                         nproc = nproc, clone = data_clone_val,
                                         split_light = split_light, ...)
  } else{
    all_germlines <- NULL
  }
  
  if(!is.na(subsample_size)){
    if(!is.numeric(subsample_size)){
      stop("subsample_size must be a numeric")
    }
    
    clones <- sampleClones(clones, size = subsample_size)
    
    if(subsample_size == 1 & !is.na(subsample_size)){
      cells <- unlist(lapply(clones$data, function(x) x@data$sequence_id))
      
      if (!is.na(cell) & cell %in% colnames(data)) {
        filtered <- data[data$sequence_id %in% cells, ]
        cells <- data$sequence_id[data[[cell]] %in% filtered[[cell]]]
      }
      
      sub_data <- data[data$sequence_id %in% cells,]
      
      clones <- formatClones(sub_data, nproc = nproc, clone = clone,
                             filterstop = TRUE, chain = chain, minseq = 1, 
                             dup_singles = T, split_light = split_light, ...)
    }
    saveRDS(clones, file = file.path(dir, "clones.rds"))
  }
  
  if(repertoire_wide){
    if(quiet > 0){
      print("constructing trees")
    }
    
    clones <- tryCatch({
      if(chain == "HL" & partition != "hl"){
        warning("Paired analysis is being requested but the paired partition is not being requested.",
                " To build the best paired trees use partition = 'hl'")
      } 
      
      getTrees(clones, build = "igphyml", exec = exec, rm_temp = FALSE, dir = dir,
               asrp = TRUE, nproc = nproc, partition = partition, ...)
    }, error = function(e){
      stop(paste0("Error during tree building:\n", 
                  e$message, "\n\n", 
                  "This may be caused by one or more problemeatic clones, \n", 
                  "Consider rerunning with `repertoire_wide = FALSE` ", 
                  "to isolate failures at the clone level.\nThis will also allow ",
                  "for results for the nonfailing clones"),
           call. = FALSE)
    })
    saveRDS(clones, file.path(dir, "clones.rds"))
  }
  
  if(quiet > 0){
    print("preparing the clones for UCA analysis")
  }
  
  clones <- invisible(do.call(rbind, parallel::mclapply(
    clones$clone_id, function(x){
    result <- tryCatch({
      processCloneGermline(clone_ids = x, clones = clones, data = data, dir = dir,
                           id = id, resolve_germ = resolve_germ, 
                           all_germlines = all_germlines, quiet = quiet, 
                           clone = clone, data_clone = data_clone_val, 
                           chain = chain, check_genes = check_genes, 
                           exec = exec, partition = partition,
                           search = search, cell = cell, 
                           split_light = split_light, ...)
    }, error = function(e){
      message(paste("Error in clone", x, ":", conditionMessage(e)))
      return(NULL)
    })
    result
  }, mc.cores = nproc)))
  
  if(nrow(clones) == 0){
    stop('No clones remain')
  }
  saveRDS(clones, file.path(dir, "clones.rds"))
  
  
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
  
  clones <- updateClone(clones = clones, data = data, references = references, 
                        dir = dir, id = id, nproc = nproc, clone_id = data_clone_val, 
                        resolve_germ = resolve_germ, fill_gaps = fill_gaps)
  saveRDS(clones, file.path(dir, "clones.rds"))
  
  unlink(rm_dir,recursive=TRUE)
  return(clones)
}
