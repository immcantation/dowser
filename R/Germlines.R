## Functions for constructing clonal germline sequences
## Based closely on CreateGermlines.py

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
createGermlines <- function(data, references, locus="locus",
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
