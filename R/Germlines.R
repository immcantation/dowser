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
findConsensus <- function(receptors, v_call = "v_call", j_call = "j_call",
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
  rec_v <- strsplit(cons[[v_call]], ",")[[1]]
  rec_v <- rec_v[rec_v %in% v_cons][1]
  rec_j <- strsplit(cons[[j_call]], ",")[[1]]
  rec_j <- rec_j[rec_j %in% j_cons][1]
  temp <- data.frame(clone_id = receptors$clone_id_unique[1], locus = receptors$locus[1], v_call = rec_v, 
                     j_call = rec_j, cons_id = cons_id)
  return(temp)
}

# find the igblast gene lengths of a given sequence 
findLengths <- function(data, sequence_id, clone, locus = "locus", seq = "sequence_id", 
                        v_len = "v_germline_length", d_len = "d_germline_length",
                        j_len = "j_germline_length", np1_len = "np1_length", 
                        np2_len = "np2_length"){
  v_length <- data[[v_len]][data[[seq]] == sequence_id]
  np1_length <- data[[np1_len]][data[[seq]] == sequence_id]
  j_length <- data[[j_len]][data[[seq]] == sequence_id]
  if(data[[locus]][data[[seq]] == sequence_id] == "IGH"){
    d_length <- data[[d_len]][data[[seq]] == sequence_id]
    np2_length <- data[[np2_len]][data[[seq]] == sequence_id]
    if(is.na(d_length)){
      d_length <- 0
    }
    if(is.na(np2_length)){
      np2_length <- 0
    }
  } else{
    d_length <- NA
    np2_length <- NA
  }
  temp <- data.frame(sequence_id = sequence_id, locus = data[[locus]][data[[seq]] == sequence_id],
                     v_length = v_length, np1_length = np1_length, d_length = d_length, 
                     np2_length = np2_length, j_length = j_length)
  # update the lengths to take out gaps -- see numbers
  numbers <- clone$data[[1]]@numbers
  regions <- clone$data[[1]]@region
  if(clone$data[[1]]@phylo_seq == "hlsequence"){
    restart_point <- which(diff(numbers) < 0) + 1
    if(data[[locus]][data[[seq]] == sequence_id] == "IGH"){
      numbers <- numbers[1:restart_point-1]
      regions <- regions[1:restart_point-1]
    }else{
      numbers <- numbers[restart_point:length(numbers)]
      regions <- regions[restart_point:length(regions)]
    }
  }
  gaps <- dplyr::setdiff(1:max(numbers), numbers)
  v_length <- 1:temp$v_length
  v_length <- length(v_length[-which(v_length %in% gaps)])
  np1_length <- (temp$v_length + 1): (temp$np1_length + temp$v_length)
  np1_length <- length(np1_length[!np1_length %in% gaps])
  if(data[[locus]][data[[seq]] == sequence_id] == "IGH"){
    d_length <- (temp$v_length + temp$np1_length + 1): (temp$d_length + temp$np1_length + temp$v_length)
    d_length <- length(d_length[!d_length %in% gaps])
    np2_length <- (temp$v_length + temp$np1_length + temp$d_length + 1): (temp$np2_length + 
                                                                            temp$d_length + temp$np1_length + temp$v_length)
    np2_length <- length(np2_length[!np2_length %in% gaps])
    j_length <- (temp$v_length + temp$np1_length + temp$d_length + temp$np2_length + 1): 
      (temp$j_length + temp$np2_length + temp$d_length + temp$np1_length + temp$v_length)
    j_length <- length(j_length[!j_length %in% gaps])
  } else{
    j_length <- (temp$v_length + temp$np1_length + 1): (temp$j_length + temp$np1_length + temp$v_length)
    j_length <- length(j_length[!j_length %in% gaps])
    d_length <- NA
    np2_length <- NA
  }
  temp$v_length <- v_length
  temp$np1_length <- np1_length
  temp$d_length <- d_length
  temp$np2_length <- np2_length
  temp$j_length <- j_length
  temp$regions <- I(list(regions))
  temp$numbers <- I(list(numbers))
  return(temp)
}

# finds where cdr3 in in a given gene
cdr3inGene <- function(temp){
  regions <- temp$regions[[1]]
  v_cdr3 <- which(regions[1:temp$v_length] == "cdr3")
  if(temp$locus == "IGH"){
    j_start <- sum(temp$v_length, temp$np1_length, temp$d_length, temp$np2_length)
    j_cdr3 <- which(regions[j_start:length(regions)] == "cdr3")
  }else{
    j_start <- sum(temp$v_length, temp$np1_length)
    j_cdr3 <- which(regions[j_start:length(regions)] == "cdr3")
  }
  temp$v_cdr3 <- I(list(v_cdr3))
  temp$j_cdr3 <- I(list(j_cdr3))
  return(temp)
}

# if check_genes is TRUE then this function goes through the proposed starting germline and 
# fixes the v and j genes of said germline 
updateGermline <- function(stats_df, uca, references, data, v_glen = "v_germline_length", 
                           np1_length = "np1_length", d_glen = "d_germline_length",
                           np2_length = "np2_length", j_glen = "j_germline_length",
                           seq = "sequence_id", germ = "germline_alignment", 
                           v_call = "v_call", j_call = "j_call", v_germline_start = "v_germline_start",
                           v_germline_end = "v_germline_end", j_germline_start = "j_germline_start",
                           j_germline_end = "j_germline_end", locus = "locus"){
  v_indx <- stats_df$v_cdr3[[1]]
  j_indx <- stats_df$j_cdr3[[1]]
  #cg_germ <- data[[germ]][data[[seq]] == stats_df[[seq]]]
  if(length(v_indx) == 0){
    v_indx <- NULL
  } else{
    if(length(v_indx) == 1 & sum(is.na(v_indx)) == 1){
      v_indx <- NULL
    }
  }
  if(length(j_indx) == 0){
    j_indx <- NULL
  } else{
    if(length(j_indx) == 1 & sum(is.na(j_indx)) == 1){
      j_indx <- NULL
    }
  }
  if(is.na(stats_df$d_length)){
    stats_df$d_length <- 0
  }
  if(is.na(stats_df$np2_length)){
    stats_df$np2_length <- 0
  }
  if(length(v_indx) > 0){
    uca <- strsplit(uca, "")[[1]]
    ref_seq <- references[[stats_df[[locus]]]]$V
    ref_seq <- ref_seq[names(ref_seq) == stats_df[[v_call]]]
    gap_values <- dplyr::setdiff(1:max(stats_df$numbers[[1]]), stats_df$numbers[[1]])
    ref_seq <- strsplit(substring(ref_seq, data[[v_germline_start]][data[[seq]] == stats_df[[seq]]],
                                  data[[v_germline_end]][data[[seq]] == stats_df[[seq]]]),
                        "")[[1]][-gap_values]
    v_indx <- v_indx[v_indx < length(ref_seq) + 1]
    uca[v_indx] <- ref_seq[v_indx]
    uca <- paste(uca, collapse = "")
  }
  if(length(j_indx) > 0){
    front_uca <- substring(uca, 1, sum(stats_df$v_length, stats_df$np1_length, stats_df$d_length, 
                                       stats_df$np2_length))
    back_uca <- substring(uca, sum(stats_df$v_length, stats_df$np1_length, stats_df$d_length, 
                                   stats_df$np2_length, 1), nchar(uca))
    back_uca <- strsplit(back_uca, "")[[1]]
    ref_seq <- references[[stats_df[[locus]]]]$J
    ref_seq <- ref_seq[names(ref_seq) == stats_df[[j_call]]]
    ref_seq <- strsplit(substring(ref_seq, data[[j_germline_start]][data[[seq]] == stats_df[[seq]]], 
                                  data[[j_germline_end]][data[[seq]] == stats_df[[seq]]]), "")[[1]]
    gap_values <- dplyr::setdiff(1:max(stats_df$numbers[[1]]), stats_df$numbers[[1]])
    front_length <- sum(data[[v_glen]][data[[seq]] == stats_df[[seq]]],
                        data[[np1_length]][data[[seq]] == stats_df[[seq]]],
                        data[[d_glen]][data[[seq]] == stats_df[[seq]]],
                        data[[np2_length]][data[[seq]] == stats_df[[seq]]])
    if(is.na(front_length)){
      front_length <- sum(data[[v_glen]][data[[seq]] == stats_df[[seq]]],
                          data[[np1_length]][data[[seq]] == stats_df[[seq]]])
    }
    gap_values <- gap_values[gap_values > front_length]
    if(length(gap_values) > 0){
      ref_seq <- ref_seq[-gap_values]
    }
    back_uca[j_indx] <- ref_seq[j_indx]
    uca <- paste0(front_uca, paste0(back_uca, collapse = ""))
    j_replacements <- j_indx + nchar(front_uca)
  } else{
    j_replacements <- j_indx
  }
  temp <- data.frame(uca = uca, v_replace = I(list(v_indx)), 
                     j_replace = I(list(j_replacements)))
  return(temp)
}

# Finds the AA of a codon
findCodon <- function(germ, group, indx, aa = FALSE){
  min_val <- group[[indx]][1]
  max_val <- group[[indx]][3]
  if(!aa){
    seq_val <- substring(germ, min_val, max_val)
  } else{
    seq_val <- alakazam::translateDNA(substring(germ, min_val, max_val)) 
  }
  return(seq_val)
}

# update the Tree Df based on the V/Js in the cdr3 boundary
updateTreeDF <- function(df, sites, codons, site = "site", codon = "codon"){
  noneffectdf <- df[!df[[site]] %in% sites,]
  temp <- df[df[[site]] %in% sites,]
  temp <- do.call(rbind, lapply(sites, function(x){
    sub <- temp[temp[[site]] == x,]
    sub <- sub[sub$codon == codons[which(sites == x)],]
    return(sub)
  }))
  df_new <- rbind(noneffectdf, temp)
  df_new <- df_new[order(df_new[[site]]),]
  return(df_new)
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
# @param resolve_v  Resolve the V gene as well?
# @param resolve_j  Resolve the J gene in addition to the V gene?
# @param all_germlines The output of createAllGermlines. Only needed if resolving genes.
# @param quiet      How much noise to print out
# @param chain      HL or H?
# @param check_genes Check if the inferred V/J lengths go into the inferred cdr3 region and adjust accordingly. 
# @param references IMGT references read in using \link{readIMGT}
#
processCloneGermline <- function(clone_ids, clones, data, dir, build, id,
                                 resolve_v = FALSE, resolve_j = FALSE,
                                 all_germlines = NULL, quiet = 0, chain = chain,
                                 check_genes = FALSE, references = NULL, ...){
  sub <- dplyr::filter(clones, !!rlang::sym("clone_id") == clone_ids)
  subDir <- file.path(dir, paste0(id, "_",clone_ids))
  if(!dir.exists(subDir)){
    dir.create(subDir, recursive = T)
  }
  saveRDS(sub, file.path(subDir, "clone.rds"))
  
  if(sub$data[[1]]@phylo_seq == "hlsequence"){
    # check to see if the lgermline and hlgermline LC are the same
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
    # need to subset regions to only the HC
    heavy_r <- r[1:nchar(sub$data[[1]]@germline)]
    light_r <- r[(nchar(sub$data[[1]]@germline) + 1): length(r)]
    cdr3_index <- (min(which(heavy_r == "cdr3")) - 3):(max(which(heavy_r == "cdr3")) + 3)
  } else if(sub$data[[1]]@phylo_seq == "lsequence"){
    cdr3_index <- (min(which(r == "cdr3")) - 3):(max(which(r == "cdr3")) + 3)
  }
  # double check that the sequence starts with C and ends with F/W
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
    # make the sample folder and save a copy of the codon table as the sample...rootprobs.txt
    # if(!dir.exists(file.path(subDir, "sample"))){
    #   dir.create(file.path(subDir, "sample"))
    # }
    file.copy(file.path(subDir, "codon_table.txt"), 
              file.path(subDir, paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")))
  }

  colnames(tree_df) = c("site", "codon", "partial_likelihood", "nope", "nada", "no", "equilbrium")
  tree_df$value <- tree_df$partial_likelihood + log(tree_df$equilbrium)
  
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
      sub_tree_df$value <- sub_tree_df$partial_likelihood + log(sub_tree_df$equilbrium)
      value <- sub_tree_df$codon[sub_tree_df$value == max(sub_tree_df$value)]
      mrcacdr3 <- paste0(value[1], substring(mrcacdr3, 4, nchar(mrcacdr3)))
    }
    if(!test_cdr3[length(test_cdr3)] %in% c("F", "W")){
      codon_site <- which(sapply(groupedList, function(group) max(cdr3_index) %in% group))
      sub_tree_df <- dplyr::filter(tree_df, !!rlang::sym("site") == codon_site - 1)
      sub_tree_df$aa <- alakazam::translateDNA(sub_tree_df$codon)
      sub_tree_df <- dplyr::filter(sub_tree_df, !!rlang::sym("aa") %in% c("W", "F"))
      sub_tree_df$value <- sub_tree_df$partial_likelihood + log(sub_tree_df$equilbrium)
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
    # get the heavy chain V and associated stats
    v_len <- min(cdr3_index)-1
    v <- substring(imgt_germline, 1, v_len)
    # get the heavy chain J and associated stats
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
  if(sub$data[[1]]@phylo_seq == "hlsequence" & resolve_v){
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
  }
  if(resolve_v){
    if(quiet > 0){
      print(paste("resolving genes for", clone_ids))
    }
    if(sub$data[[1]]@phylo_seq == "sequence"){
      has_multiple <- all_germlines[all_germlines$clone_id == clone_ids,]
      heavy_indx <- grepl("^IGH", has_multiple$v_call)
      has_multiple <- has_multiple[heavy_indx,]
    }
    # add the padding to the ungapped
    has_multiple$ungapped <- unlist(parallel::mclapply(1:nrow(has_multiple), function(y){
      current <- has_multiple$ungapped[y]
      if(nchar(current) %% 3 > 0){
        current <- paste0(current, paste0(rep("N", 3 - nchar(current) %% 3), collapse = ""))
      }
      return(current)
    }))
    saveRDS(has_multiple, file.path(subDir, "all_germlines_heavy.rds"))
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
                                !!rlang::sym("site") %in% c((v_stop + nchar(mrcacdr3)/3): max(tree_df$site)))
      vj <- paste0(v_alt, j_alt, collapse = "")
      gene_list <- strsplit(vj, "")[[1]]
      groupedSeq <- split(gene_list, ceiling(seq_along(gene_list) / 3))
      if("N" %in% groupedSeq[[length(groupedSeq)]]){
        # add the N option to the J
        df_row <- sub_df[sub_df$site == max(sub_df$site),]
        df_row$value <- df_row$partial_likelihood + log(df_row$equilbrium)
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
    if(resolve_j){
      j <- germlines$j
    }
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
    
    # fill the nps with mrca but make sure conserved sites are still there 
    light_cdr3_test <- strsplit(alakazam::translateDNA(light_cdr3), "")[[1]]
    
    if(light_cdr3_test[1] != "C" || !light_cdr3_test[length(light_cdr3_test)] %in% c("F", "W")){
      if(light_cdr3_test[1] != "C"){
        codon_site <- nchar(v_light)/3 + 1
        sub_tree_df <- dplyr::filter(tree_df_light, !!rlang::sym("site") == codon_site)
        sub_tree_df$aa <- alakazam::translateDNA(sub_tree_df$codon)
        sub_tree_df <- dplyr::filter(sub_tree_df, !!rlang::sym("aa") == "C")
        sub_tree_df$value <- sub_tree_df$partial_likelihood + log(sub_tree_df$equilbrium)
        value <- sub_tree_df$codon[sub_tree_df$value == max(sub_tree_df$value)]
        light_cdr3 <- paste0(value[1], substring(light_cdr3, 4, nchar(light_cdr3)))
      }
      if(!light_cdr3_test[length(light_cdr3_test)] %in% c("F", "W")){
        codon_site <- nchar(paste0(v_light, light_cdr3))/3 + 1
        sub_tree_df <- dplyr::filter(tree_df_light, !!rlang::sym("site") == codon_site)
        sub_tree_df$aa <- alakazam::translateDNA(sub_tree_df$codon)
        sub_tree_df <- dplyr::filter(sub_tree_df, !!rlang::sym("aa") %in% c("W", "F"))
        sub_tree_df$value <- sub_tree_df$partial_likelihood + log(sub_tree_df$equilbrium)
        value <- sub_tree_df$codon[sub_tree_df$value == max(sub_tree_df$value)]
        light_cdr3 <- paste0(substring(light_cdr3, 1, nchar(light_cdr3)-3), value[1])
      }
    }
    if(resolve_v){
      germlines_light <- do.call(rbind, parallel::mclapply(1:nrow(has_multiple_light), function(x){
        value <- has_multiple_light$ungapped[x]
        v_alt <- substring(value, 1, nchar(v_light))
        j_alt <- substring(value, nchar(paste0(v_light, light_cdr3)) + 1, nchar(value))
        if(nchar(j_alt) > nchar(j_light)){
          j_alt <- substring(j_alt, 1, nchar(j_light))
        }
        sub_tree_df <- tree_df_light[tree_df_light$site %in% c(1:(nchar(v_light)/3),
                        (max(tree_df_light$new_site) - (nchar(j_light)/3) + 1):max(tree_df_light$new_site)),]
        # get the likelihood of the v_alt and j_alt combo
        vj <- paste0(v_alt, j_alt, collapse = "")
        gene_list <- strsplit(vj, "")[[1]]
        groupedSeq <- split(gene_list, ceiling(seq_along(gene_list) / 3))
        if("N" %in% groupedSeq[[length(groupedSeq)]]){
          # add the N option to the J
          df_row <- sub_tree_df[sub_tree_df$site == max(sub_tree_df$site),]
          df_row$value <- df_row$partial_likelihood + log(df_row$equilbrium)
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
      if(resolve_j){
        j_light <- germlines_light$j
      }
    } 
  }
  saveRDS(sub, file.path(subDir, "clone.rds"))
  # check if the V/J goes into the CDR3 region definition 
  if(check_genes){
    if(sub$data[[1]]@phylo_seq == "hlsequence"){
      # use heavy and light_r
      heavy_cons <- findConsensus(dplyr::filter(data, !!rlang::sym("clone_id") == sub$clone_id &
                                                  !!rlang::sym('locus') == "IGH"))
      light_cons <- findConsensus(dplyr::filter(data, !!rlang::sym("clone_id") == sub$clone_id &
                                                  !!rlang::sym('locus') != "IGH"))
      heavy_stats <- findLengths(data, heavy_cons$cons_id, sub)
      light_stats <- findLengths(data, light_cons$cons_id, sub)
      heavy_stats <- cdr3inGene(heavy_stats)
      light_stats <- cdr3inGene(light_stats)
      # grab the v and j calls 
      heavy_stats$v_call <- heavy_cons$v_call
      heavy_stats$j_call <- heavy_cons$j_call
      light_stats$v_call <- light_cons$v_call
      light_stats$j_call <- light_cons$j_call
      
      heavy_germ <- updateGermline(heavy_stats, paste0(v, mrcacdr3, j), references,
                                  data)
      light_germ <- updateGermline(light_stats, paste0(v_light, light_cdr3, j_light),
                                   references, data)
      # now update the tree df to only look for the codon combos for the replaced sites 
      # are AAs of what is there 
      heavy_groups <- lapply(seq(1, nchar(sub$data[[1]]@germline), by = 3),
                             function(i) i:min(i+2, nchar(sub$data[[1]]@germline)))
      light_groups <- lapply(seq(1, nchar(sub$data[[1]]@lgermline), by = 3),
                             function(i) i:min(i+2, nchar(sub$data[[1]]@lgermline)))

      if(!is.null(heavy_germ$v_replace[[1]]) | !is.null(heavy_germ$j_replace[[1]])){
        # find the sites that were changed
        heavy_indices <- unlist(unique(lapply(c(heavy_germ$v_replace[[1]], heavy_germ$j_replace[[1]]), function(val) {
          which(sapply(heavy_groups, function(g) val %in% g))})))
        # find the AA of what should be there 
        heavy_codon <- unlist(lapply(1:length(heavy_indices), function(x){
          value <- findCodon(heavy_germ$uca, heavy_groups, heavy_indices[x])
          return(value)
        }))
        # if there are Xs or *s in the returned values remove those indices 
        stop_codons <- c("TAA", "TAG", "TGA")
        heavy_indices <- heavy_indices[!heavy_codon %in% stop_codons]
        heavy_codon <- heavy_codon[!heavy_codon %in% stop_codons]
        heavy_codon <- heavy_codon[!heavy_indices %in% c(nchar(v)/3, nchar(paste0(v, mrcacdr3))/3)]
        heavy_indices <- heavy_indices[!heavy_indices %in% c(nchar(v)/3, nchar(paste0(v, mrcacdr3))/3)]
        # remove any indices that have codons that translate to X
        heavy_aa <- unlist(lapply(1:length(heavy_indices), function(x){
          value <- findCodon(heavy_germ$uca, heavy_groups, heavy_indices[x], aa= TRUE)
          return(value)
        }))
        heavy_indices <- heavy_indices[heavy_aa != "X"]
        heavy_codon[heavy_aa != "X"]
        # update to be 0 based counting 
        heavy_indices <- heavy_indices - 1
        # filter down the associated df 
        tree_df <- updateTreeDF(tree_df, heavy_indices, heavy_codon)
        # if there is a stop codon find it and replace with what the og value was
        germ_test <- which(strsplit(alakazam::translateDNA(heavy_germ$uca), "")[[1]] == "*")
        if(length(germ_test) > 0){
          change_vals <- heavy_groups[[germ_test]]
          uca <- strsplit(heavy_germ$uca, "")[[1]]
          og_uca <- strsplit(paste0(v, mrcacdr3, j), "")[[1]]
          uca[change_vals] <- og_uca[change_vals]
          uca <- paste(uca, collapse = "")
          heavy_germ$uca <- uca
        }
        # update the v, cdr3, and j cuts 
        v <- substring(heavy_germ$uca, 1, nchar(v))
        cdr3_c <- substring(mrcacdr3, 1, 3)
        cdr3_fw <- substring(mrcacdr3, nchar(mrcacdr3)-2, nchar(mrcacdr3))
        mrcacdr3 <- paste0(cdr3_c, substring(heavy_germ$uca, nchar(v) + 4, nchar(v) + nchar(mrcacdr3)-3),
                           cdr3_fw)
        j <- substring(heavy_germ$uca, nchar(v) + nchar(mrcacdr3) + 1, nchar(heavy_germ$uca))
      }
      if(!is.null(light_germ$v_replace[[1]]) | !is.null(light_germ$j_replace[[1]])){
        light_indices <- unlist(unique(lapply(c(light_germ$v_replace[[1]], light_germ$j_replace[[1]]), function(val) {
          which(sapply(light_groups, function(g) val %in% g))})))
        light_codon <- unlist(lapply(1:length(light_indices), function(x){
          value <- findAA(light_germ$uca, light_groups, light_indices[x])
          return(value)
        }))
        stop_codons <- c("TAA", "TAG", "TGA")
        light_indices <- light_indices[!light_codon %in% stop_codons]
        light_codon <- light_codon[!lilight_codonght_aa %in% stop_codons]
        light_codon <- light_codon[!light_indices %in% c(nchar(v_light)/3, nchar(paste0(v_light, light_cdr3))/3)]
        light_indices <- light_indices[!light_indices %in% c(nchar(v_light)/3, nchar(paste0(v_light, light_cdr3))/3)]
        light_aa <- unlist(lapply(1:length(light_indices), function(x){
          value <- findCodon(light_germ$uca, light_groups, light_indices[x], aa= TRUE)
          return(value)
        }))
        light_indices <- light_indices[light_aa != "X"]
        light_codon[light_aa != "X"]
        light_indices <- light_indices - 1
        tree_df_light <- updateTreeDF(tree_df_light, light_indices, light_codon)
        germ_test <- which(strsplit(alakazam::translateDNA(light_germ$uca), "")[[1]] == "*")
        if(length(germ_test) > 0){
          change_vals <- light_groups[[germ_test]]
          uca <- strsplit(light_germ$uca, "")[[1]]
          og_uca <- strsplit(paste0(v_light, light_cdr3, j_light), "")[[1]]
          uca[change_vals] <- og_uca[change_vals]
          uca <- paste(uca, collapse = "")
          light_germ$uca <- uca
        }
        v_light <- substring(light_germ$uca, 1, nchar(v_light))
        cdr3_c <- substring(light_cdr3, 1, 3)
        cdr3_fw <- substring(light_cdr3, nchar(light_cdr3)-2, nchar(light_cdr3))
        light_cdr3 <- paste0(cdr3_c, substring(light_germ$uca, nchar(v_light) + 4, nchar(v_light) + nchar(light_cdr3)-3),
                             cdr3_fw)
        j_light <- substring(light_germ$uca, nchar(v_light) + nchar(light_cdr3) + 1, nchar(light_germ$uca))
      }
    }else if(sub$data[[1]]@phylo_seq == "sequence"){
      heavy_cons <- findConsensus(dplyr::filter(data, !!rlang::sym("clone_id") == sub$clone_id &
                                                  !!rlang::sym('locus') == "IGH"))
      heavy_stats <- findLengths(data, heavy_cons$cons_id, sub)
      heavy_stats <- cdr3inGene(heavy_stats)
      heavy_stats$v_call <- heavy_cons$v_call
      heavy_stats$j_call <- heavy_cons$j_call
      heavy_germ <- updateGermline(heavy_stats, paste0(v, mrcacdr3, j), references,
                                   data)
      if(!is.null(heavy_germ$v_replace[[1]]) | !is.null(heavy_germ$j_replace[[1]])){
        heavy_groups <- lapply(seq(1, nchar(sub$data[[1]]@germline), by = 3),
                               function(i) i:min(i+2, nchar(sub$data[[1]]@germline)))
        heavy_indices <- unlist(unique(lapply(c(heavy_germ$v_replace[[1]], heavy_germ$j_replace[[1]]), function(val) {
          which(sapply(heavy_groups, function(g) val %in% g))})))
        heavy_codon <- unlist(lapply(1:length(heavy_indices), function(x){
          value <- findCodon(heavy_germ$uca, heavy_groups, heavy_indices[x])
          return(value)
        }))
        heavy_indices <- heavy_indices[!heavy_codon %in% stop_codons]
        heavy_codon <- heavy_codon[!heavy_codon %in% stop_codons]
        heavy_codon <- heavy_codon[!heavy_indices %in% c(nchar(v)/3, nchar(paste0(v, mrcacdr3))/3)]
        heavy_indices <- heavy_indices[!heavy_indices %in% c(nchar(v)/3, nchar(paste0(v, mrcacdr3))/3)]
        heavy_aa <- unlist(lapply(1:length(heavy_indices), function(x){
          value <- findCodon(heavy_germ$uca, heavy_groups, heavy_indices[x], aa= TRUE)
          return(value)
        }))
        heavy_indices <- heavy_indices[heavy_aa != "X"]
        heavy_codon[heavy_aa != "X"]
        heavy_indices <- heavy_indices - 1
        tree_df <- updateTreeDF(tree_df, heavy_indices, heavy_codon)
        germ_test <- which(strsplit(alakazam::translateDNA(heavy_germ$uca), "")[[1]] == "*")
        if(length(germ_test) > 0){
          change_vals <- heavy_groups[[germ_test]]
          uca <- strsplit(heavy_germ$uca, "")[[1]]
          og_uca <- strsplit(paste0(v, mrcacdr3, j), "")[[1]]
          uca[change_vals] <- og_uca[change_vals]
          uca <- paste(uca, collapse = "")
          heavy_germ$uca <- uca
        }
        v <- substring(heavy_germ$uca, 1, nchar(v))
        cdr3_c <- substring(mrcacdr3, 1, 3)
        cdr3_fw <- substring(mrcacdr3, nchar(mrcacdr3)-2, nchar(mrcacdr3))
        mrcacdr3 <- paste0(cdr3_c, substring(heavy_germ$uca, nchar(v) + 4, nchar(v) + nchar(mrcacdr3)-3),
                           cdr3_fw)
        j <- substring(heavy_germ$uca, nchar(v) + nchar(mrcacdr3) + 1, nchar(heavy_germ$uca))
        # save tree_df
        tree_df <- tree_df[, !names(tree_df) == "value"]
        write.table(tree_df, file.path(subDir, paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")), 
                    quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
      }
      
    } else if(sub$data[[1]]@phylo_seq == "lsequence"){
      light_cons <- findConsensus(dplyr::filter(data, !!rlang::sym("clone_id") == sub$clone_id &
                                                  !!rlang::sym('locus') != "IGH"))
      light_stats <- findLengths(data, light_cons$cons_id, sub)
      light_stats <- cdr3inGene(light_stats)
      light_stats$v_call <- light_cons$v_call
      light_stats$j_call <- light_cons$j_call
      light_germ <- updateGermline(light_stats, paste0(v, mrcacdr3, j),
                                   references, data)
      if(!is.null(light_germ$v_replace[[1]]) | !is.null(light_germ$j_replace[[1]])){
        light_groups <- lapply(seq(1, nchar(sub$data[[1]]@lgermline), by = 3),
                               function(i) i:min(i+2, nchar(sub$data[[1]]@lgermline)))
        light_indices <- unlist(unique(lapply(c(light_germ$v_replace[[1]], light_germ$j_replace[[1]]), function(val) {
          which(sapply(light_groups, function(g) val %in% g))})))
        light_codon <- unlist(lapply(1:length(light_indices), function(x){
          value <- findAA(light_germ$uca, light_groups, light_indices[x])
          return(value)
        }))
        stop_codons <- c("TAA", "TAG", "TGA")
        light_indices <- light_indices[!light_codon %in% stop_codons]
        light_codon <- light_codon[!lilight_codonght_aa %in% stop_codons]
        light_aa <- light_aa[!light_indices %in% c(nchar(v)/3, nchar(paste0(v, mrcacdr3))/3)]
        light_indices <- light_indices[!light_indices %in% c(nchar(v)/3, nchar(paste0(v, mrcacdr3))/3)]
        light_aa <- unlist(lapply(1:length(light_indices), function(x){
          value <- findCodon(light_germ$uca, light_groups, light_indices[x], aa= TRUE)
          return(value)
        }))
        light_indices <- light_indices[light_aa != "X"]
        light_codon[light_aa != "X"]
        light_indices <- light_indices - 1
        tree_df <- updateTreeDF(tree_df, light_indices, light_codon)
        germ_test <- which(strsplit(alakazam::translateDNA(light_germ$uca), "")[[1]] == "*")
        if(length(germ_test) > 0){
          change_vals <- light_groups[[germ_test]]
          uca <- strsplit(light_germ$uca, "")[[1]]
          og_uca <- strsplit(paste0(v, mrcacdr3, j), "")[[1]]
          uca[change_vals] <- og_uca[change_vals]
          uca <- paste(uca, collapse = "")
          light_germ$uca <- uca
        }
        v <- substring(light_germ$uca, 1, nchar(v))
        cdr3_c <- substring(mrcacdr3, 1, 3)
        cdr3_fw <- substring(mrcacdr3, nchar(mrcacdr3)-2, nchar(mrcacdr3))
        mrcacdr3 <- paste0(cdr3_c, substring(light_germ$uca, nchar(v) + 4, nchar(v) + nchar(mrcacdr3)-3),
                           cdr3_fw)
        j <- substring(light_germ$uca, nchar(v) + nchar(mrcacdr3) + 1, nchar(light_germ$uca))
        tree_df <- tree_df[, !names(tree_df) == "value"]
        write.table(tree_df, file.path(subDir, paste0(clone_ids, ".fasta_igphyml_rootprobs_hlp.txt")), 
                    quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
      }
      
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
    # if(!dir.exists(file.path(subDir, "sample"))){
    #   dir.create(file.path(subDir, "sample"))
    # }
    tree_df <- tree_df[, !names(tree_df) %in% c("value")]
    write.table(tree_df, file.path(subDir, "heavy_table.txt"), quote = FALSE,
                sep = "\t", col.names = FALSE, row.names = FALSE)
    tree_df_light <- tree_df_light[, !names(tree_df_light) %in% c("value")]
    write.table(tree_df_light, file.path(subDir, "light_table.txt"), quote = FALSE,
                sep = "\t", col.names = FALSE, row.names = FALSE)
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
# @param method     The method to use with OLGA. Options are ML (maximum likelihood) or MCMC
# @param tree_specs A string that contains the igphyml exec, rate, motif, hotness, and omega values used for tree building separated by ;
# @param rscript    The file path to rsrcipt needed to rebuild trees for ML_tree. 
# @param clone      The column name for the clone identifier 
#

callOlga <- function(clones, dir, uca_script, python, max_iters, nproc, id, model_folder,
                     model_folder_igk, model_folder_igl, quiet, method, tree_specs, 
                     rscript = NULL, clone = "clone_id"){
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
    "--directory", dir, 
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
    "--method", method, 
    "--tree_specs", tree_specs,
    "--rscript", ifelse(is.null(rscript), "NULL", path.expand(rscript))
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
      buildAllClonalGermlines(
        receptors = sub[sub[[locus]] == l,],
        references = references,
        chain = l, 
        seq = seq, 
        v_call = v_call, 
        d_call = d_call,
        j_call = j_call, 
        amino_acid = amino_acid,
        id = id, 
        clone = clone, 
        v_germ_start = v_germ_start,
        v_germ_end = v_germ_end,
        v_germ_length = v_germ_length,
        d_germ_start = d_germ_start, 
        d_germ_end = d_germ_end, 
        d_germ_length = d_germ_length, 
        j_germ_start = j_germ_start, 
        j_germ_end = j_germ_end, 
        j_germ_length = j_germ_length,
        np1_length = np1_length, 
        np2_length = np2_length, ...)
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
    }
    masked_sites <- which(comp_df$diff)
    clone_germ[masked_sites]  <- "N"
    if(chain == "H"){
      sub$data[[1]]@germline <- paste0(clone_germ, collapse = "")
    } else if(chain == "HL"){
      sub$data[[1]]@hlgermline <- paste0(clone_germ, collapse = "")
    }
    return(sub)
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
#' @param omega         Omega parameters to estimate (see IgPhyMl docs)
#' @param optimize      Optimize HLP rates (r), lengths (l), and/or topology (r)
#' @param motifs        Motifs to consider (see IgPhyMl docs)
#' @param hotness       Hotness parameters to estimate (see IgPhyML docs)
#' @param resolve_v     Resolve the V gene as well 
#' @param resolve_j     Resolve the J gene as well  
#' @param references    Reference genes. See \link{readIMGT}
#' @param clone         The name of the clone id column used in \link{formatClones}
#' @param heavy         The name of the heavy chain locus. Default is IGH. 
#' @param sampling_method How to subsample. Methods include 'random', 'lm' (least mutated),
#'                      and 'ratio' (a weighted sampling with heavier weights
#'                      towards the least mutated sequences). The later two methods
#'                      require 'mu_freq' to be passed as a trait when running
#'                      \link{formatClones}
#' @param subsample_size The amount that the clone should be sampled down to. Default is 5. 
#' @param method        The method to find the UCA with. Options are ML (maximum likelihood), ML_tree, or MCMC.
#' @param rscript       The file path to the rscript with get_UCA.py
#' @param check_genes   Check if the inferred V/J lengths go into the inferred cdr3 region and adjust accordingly. 
#' @param ...           Additional arguments passed to \link{buildGermline}
#' @return An \code{airrClone} object with trees and the inferred UCA
#' @details Return object adds/edits following columns:
#' \itemize{
#'   \item  \code{trees}:  The phylogenies associated with each clone
#'   \item  \code{UCA}:    The inferred UCA
#' }
#' @seealso \link{getTrees} 
#' @export
getTreesAndUCAs <- function(clones, data, dir = NULL, build, exec,  model_folder,
                           model_folder_igk = NULL, model_folder_igl = NULL, uca_script,
                           python = "python", id = "sample", max_iters = 100, nproc = 1,
                           rm_temp = TRUE, quiet = 0, chain = "H", omega = NULL, optimize = "lr",
                           motifs = "FCH", hotness = "e,e,e,e,e,e", resolve_v = FALSE,
                           resolve_j = FALSE, references = NULL, clone = "clone_id", heavy = "IGH",
                           sampling_method = c('random', 'lm', 'ratio'), subsample_size = 5, 
                           method = "ML", rscript = NULL, check_genes = TRUE,...){
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
    stop("check_genes cannot be run without references. Pass a references object using references =",
         "References need to be read in using dowser::readIMGT()")
  }
  
  # make the tree_specs value 
  tree_specs <- paste0(
    if (is.null(exec)) "NULL" else exec, ";",
    if (is.null(optimize)) "NULL" else optimize, ";",
    if (is.null(motifs)) "NULL" else motifs, ";",
    if (is.null(hotness)) "NULL" else hotness, ";",
    if (is.null(omega)) "NULL" else omega
  )
  
  # subsample the clones 
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
  
  if(resolve_v | resolve_j){
    if(is.null(data)){
      stop("The data object is required with current settings")
    }
    if(is.null(references)){
      stop("References must be supplied with current settings")
    }
    all_germlines <- createAllGermlines(data = data, references = references,
                                        nproc = nproc, clone = clone, trim_lengths = TRUE)
    saveRDS(all_germlines, file.path(dir, "all_germlines.rds"))
    if(resolve_v | resolve_j){
      clones <- maskAmbigousReferenceSites(clones = clones, all_germlines = all_germlines,
                                           nproc = nproc, clone = clone, chain = chain)
    }
  }
  
  if(quiet > 0){
    print("constructing trees")
  }
  if(build == "igphyml"){
    clones <- getTrees(clones, build = build, exec = exec, rm_temp = FALSE, dir = dir,
                    omega = omega, optimize = optimize, motifs = motifs, 
                    hotness = hotness, asrp = TRUE, chain = chain, ...)
  } else if(build == "pml"){
    clones <- getTrees(clones, build = build, rm_temp = FALSE, dir = dir, asrp = TRUE, ...)
  } else{
    stop("the tree bulding method ", build, "is not supported")
  }
  if(quiet > 0){
    print("preparing the clones for UCA analysis")
  }
  saveRDS(clones, file.path(dir, "clones.rds"))
  clones <- do.call(rbind, invisible(parallel::mclapply(unique(clones[[clone]]), function(x)
    processCloneGermline(clone_ids = x, clones = clones, data = data, dir = dir, build = build, 
                         id = id, resolve_v = resolve_v, resolve_j = resolve_j,
                         all_germlines = all_germlines, quiet = quiet, chain = chain, 
                         check_genes = check_genes, references = references), mc.cores = nproc)))
  # run the UCA
  if(quiet > 0){
    print("running UCA analysis")
  }
  #clone_ids_str <- paste(unique(clones$clone_id), collapse = ",")
  callOlga(clones = clones, dir = dir, model_folder = model_folder,
           model_folder_igk = model_folder_igk, model_folder_igl = model_folder_igl,
           uca_script = uca_script, python = python, max_iters = max_iters, nproc = nproc,
           id = id, quiet = quiet, method = method, tree_specs = tree_specs, rscript = rscript, 
           clone = clone)
  
  # read in the clones to make the base clones object again with the UCA in the data table
  if(quiet > 0){
    print("updating clones")
  }
  clones <- updateClone(clones, dir, id, nproc)
  
  # unlink if desired 
  unlink(rm_dir,recursive=TRUE)
  return(clones)
}
