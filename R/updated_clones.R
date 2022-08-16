# This is the test script to update the makeAirrClone function in dowser. 
# the data used is from the flu data that Ken has access to so I can ensure that we have light chains

#library(alakazam)
#library(dplyr)
# read in a change-o input and make it light chain 
#heavy <- readChangeoDb("~/Documents/Flu Data/heavy.tsv") 
#light <- readChangeoDb("~/Documents/Flu Data/light.tsv") 
# make sure there aren't double heavy chains
#heavy_count = table(heavy$cell_id)
#sum(heavy_count > 1)
#multi_heavy_cells = names(heavy_count[heavy_count > 1])
#heavy = filter(heavy, !cell_id %in% multi_heavy_cells)

# get the subclones
#comb = dowser::getSubclones(heavy,light,cell="cell_id", nproc= 5)
#comb = comb[which(!is.na(comb$c_call)),]
# since I already made this bad boy I'm going to down sample this to make it go quicker
#to_keep <- sample(comb$cell_id, 500)
#comb_sub <- filter(comb, cell_id %in% to_keep)

# get the germlines
#references = dowser::readIMGT(dir = "~/share/germlines/imgt/human/vdj")
#h = dowser::createGermlines(filter(comb_sub,locus=="IGH"),references, clone = "clone_id", nproc = 5)
#k = dowser::createGermlines(filter(comb_sub,locus=="IGK"),references,locus="IGK", clone = "clone_id", nproc = 5)
#l = dowser::createGermlines(filter(comb_sub,locus=="IGL"),references,locus="IGL", clone = "clone_id", nproc = 5)
#comb_germline = bind_rows(h, l, k)
#comb_germline = comb_germline[which(!is.na(comb_germline$c_call)),]


# set the class so this doesn't break when testing
setClass("airrClone",
         slots=c(
           data="data.frame",
           clone="character",
           germline="character", 
           lgermline="character",
           hlgermline="character",
           v_gene="character", 
           j_gene="character", 
           junc_len="numeric",
           locus="character",
           #chain="character",
           region="character",
           phylo_seq="character",
           numbers="numeric"))

# now this will eventually loop but for now I'm just trying to get this to work once so all the options are laid out here. 
data <- filter(comb_germline, clone_id == "42506")
#data <- filter(data, locus == "IGH")

# if it wasn't clear, the "L" option needs to have NO heavy chains otherwise I tell it to gets transfered to the "HL" option

# Clone processing functions
#' Generate a airrClone object for lineage construction
#' 
#' \code{makeAirrClone} takes a data.frame with AIRR or Change-O style columns as input and 
#' masks gap positions, masks ragged ends, removes duplicates sequences, and merges 
#' annotations associated with duplicate sequences. It returns a \code{airrClone} 
#' object which serves as input for lineage reconstruction.
#' 
#' @param    data         data.frame containing the AIRR or Change-O data for a clone. See Details
#'                        for the list of required columns and their default values.
#' @param    id           name of the column containing sequence identifiers.
#' @param    seq          name of the column containing observed DNA sequences. All 
#'                        sequences in this column must be multiple aligned.
#' @param    germ         name of the column containing germline DNA sequences. All entries 
#'                        in this column should be identical for any given clone, and they
#'                        must be multiple aligned with the data in the \code{seq} column.
#' @param    v_call        name of the column containing V-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    j_call        name of the column containing J-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    junc_len     name of the column containing the length of the junction as a 
#'                        numeric value. All entries in this column should be identical 
#'                        for any given clone.
#' @param    clone        name of the column containing the identifier for the clone. All 
#'                        entries in this column should be identical.
#' @param    mask_char    character to use for masking and padding.
#' @param    max_mask     maximum number of characters to mask at the leading and trailing
#'                        sequence ends. If \code{NULL} then the upper masking bound will 
#'                        be automatically determined from the maximum number of observed 
#'                        leading or trailing Ns amongst all sequences. If set to \code{0} 
#'                        (default) then masking will not be performed.
#' @param    pad_end      if \code{TRUE} pad the end of each sequence with \code{mask_char}
#'                        to make every sequence the same length.
#' @param    text_fields  text annotation columns to retain and merge during duplicate removal.
#' @param    num_fields   numeric annotation columns to retain and sum during duplicate removal.
#' @param    seq_fields   sequence annotation columns to retain and collapse during duplicate 
#'                        removal. Note, this is distinct from the \code{seq} and \code{germ} 
#'                        arguments, which contain the primary sequence data for the clone
#'                        and should not be repeated in this argument.
#' @param    add_count    if \code{TRUE} add an additional annotation column called 
#'                        \code{COLLAPSE_COUNT} during duplicate removal that indicates the 
#'                        number of sequences that were collapsed.
#' @param    verbose      passed on to \code{collapseDuplicates}. If \code{TRUE}, report the 
#'                        numbers of input, discarded and output sequences; otherwise, process
#'                        sequences silently.                        
#' @param    collapse     collapse identical sequences?
#' @param    traits       column ids to keep distinct during sequence collapse 
#' @param    chain        if HL, include light chain information if available. If L, only uses 
#'                        light chain information. 
#' @param    heavy        name of heavy chain locus (default = "IGH")
#' @param    cell         name of the column containing cell assignment information
#' @param    locus        name of the column containing locus information
#' @param    mod3         pad sequences to length mutliple three?
#' @param    randomize    randomize sequence order? Important if using PHYLIP
#' @param    use_regions   assign CDR/FWR regions?
#' @param    dup_singles   Duplicate sequences in singleton clones to include them as trees?
#' @return   A \link{airrClone} object containing the modified clone.
#'
#' @details
#' The input data.frame (\code{data}) must columns for each of the required column name 
#' arguments: \code{id}, \code{seq}, \code{germ}, \code{v_call}, \code{j_call}, 
#' \code{junc_len}, and \code{clone}.  

#' Additional annotation columns specified in the \code{traits}, \code{text_fields}, 
#' \code{num_fields} or \code{seq_fields} arguments will be retained in the \code{data} 
#' slot of the return object, but are not required. These options differ by their behavior
#' among collapsed sequences. Identical sequences that differ by any values specified in the
#' \code{traits} option will be kept distinct. Identical sequences that differ only by
#' values in the \code{num_fields} option will be collapsed and the values of their 
#' \code{num_fields} columns will be added together. Similar behavior occurs with 
#' \code{text_fields} but the unique values will concatenated with a comma.
#' 
#' The default columns are IMGT-gapped sequence columns, but this is not a requirement. 
#' However, all sequences (both observed and germline) must be multiple aligned using
#' some scheme for both proper duplicate removal and lineage reconstruction. 
#'
#' The value for the germline sequence, V-segment gene call, J-segment gene call, 
#' junction length, and clone identifier are determined from the first entry in the 
#' \code{germ}, \code{v_call}, \code{j_call}, \code{junc_len} and \code{clone} columns, 
#' respectively. For any given clone, each value in these columns should be identical.
#'
#' To allow for cases where heavy and light chains are used, this function returns three
#' sequence columns for heavy chains (sequence), light chain (lsequence, empty if none 
#' available), and concatenated heavy+light chain (hlsequence). These contain sequences
#' in alignment with germline, lgermline, and hlgermline slots, respectively. The sequence
#' column used for build trees is specified in the \code{phylo_seq} slot. Importantly, 
#' this column is also the sequence column that also has uninformative columns removed
#' by \code{cleanAlignment}. It is highly likely we will change this system to a single 
#' \code{sequence} and \code{germline} slot in the near future.
#'
#' The airrClone object also contains vectors \code{locus}, \code{region}, and 
#' \code{numbers}, which contain the locus, IMGT region, and IMGT number for each position
#' in the sequence column specified in \code{phylo_seq}. If IMGT-gapped sequences are not 
#' supplied, this will likely result in an error. Specify \code{use_regions=FALSE} if not
#' using IMGT-gapped sequences
#'  
#' @seealso  Returns an \link{airrClone}. See \link{formatClones} to generate an 
#' ordered list of airrClone objects.
#' @examples
#' data(ExampleAirr)
#' airr_clone <- makeAirrClone(ExampleAirr[ExampleAirr$clone_id=="3184",])
#' @export

makeAirrClone <- 
  function(data, id="sequence_id", seq="sequence_alignment", 
           germ="germline_alignment_d_mask", v_call="v_call", j_call="j_call",
           junc_len="junction_length", clone="clone_id", mask_char="N",
           max_mask=0, pad_end=TRUE, text_fields=NULL, num_fields=NULL, seq_fields=NULL,
           add_count=TRUE, verbose=FALSE, collapse=TRUE, chain="H", heavy=NULL,
           cell="cell_id", locus="locus", traits=NULL, mod3=TRUE, randomize=TRUE,
           use_regions=TRUE, dup_singles=FALSE){

  check <- alakazam::checkColumns(data, 
                                  unique(c(id, seq, germ, v_call, j_call, junc_len, clone, 
                                           text_fields, num_fields, seq_fields, traits)))
  if(chain=="HL"){
    check <- alakazam::checkColumns(data, c(cell,locus))
    if (check != TRUE) { stop(check) }
    if(is.null(heavy)){
      stop(paste("clone",unique(dplyr::pull(data,clone)),
                 "heavy chain loci ID must be specified if combining loci!"))
    }
    
    heavycount = max(table(data[data[[locus]] == heavy,][[cell]]))
    if(max(heavycount) > 1){
      stop(paste0(sum(heavycount > 1),
                  " cells with multiple heavy chains found. Remove before proceeeding"))
    }
    
    # Ensure cell and loci columns are not duplicated
    text_fields <- text_fields[text_fields != rlang::sym(cell)]
    text_fields <- text_fields[text_fields != rlang::sym(locus)]
    seq_fields <- seq_fields[seq_fields != rlang::sym(cell)]
    seq_fields <- seq_fields[seq_fields != rlang::sym(locus)]
    # Replace gaps with Ns and masked ragged ends
    tmp_df <- data[, unique(c(id, seq, junc_len, text_fields, num_fields, 
                              seq_fields, cell, locus, traits))]
    tmp_df[[seq]] <- alakazam::maskSeqGaps(tmp_df[[seq]], mask_char=mask_char, 
                                           outer_only=FALSE)
    hc <- dplyr::filter(tmp_df,!!rlang::sym(locus)==rlang::sym(heavy))
    alt <- dplyr::filter(tmp_df,!!rlang::sym(locus)!=rlang::sym(heavy))
    if(nrow(hc) == 0){
      stop(paste("clone",unique(dplyr::pull(data,clone)),
                 "heavy chain locus not found in dataset!"))
    }
    if(nrow(alt) == 0){
      chain <- "H"
    }else{
      if(length(unique(dplyr::pull(alt,!!locus))) > 1){
        stop(paste("clone",paste(unique(dplyr::pull(data,clone)),collapse=""),
                   "currently only one alternate loci per clone supported"))
      }
    }
  }
  if(chain=="L"){
    check <- alakazam::checkColumns(data, c(cell,locus))
    if (check != TRUE) { stop(check) }
  
  #  heavycount = max(table(data[data[[locus]] == heavy,][[cell]]))
   # if(max(heavycount) > 1){
    #  stop(paste0(sum(heavycount > 1),
     #             " cells with multiple heavy chains found. Remove before proceeeding"))
    #}
    
    # Ensure cell and loci columns are not duplicated
    text_fields <- text_fields[text_fields != rlang::sym(cell)]
    text_fields <- text_fields[text_fields != rlang::sym(locus)]
    seq_fields <- seq_fields[seq_fields != rlang::sym(cell)]
    seq_fields <- seq_fields[seq_fields != rlang::sym(locus)]
    # Replace gaps with Ns and masked ragged ends
    tmp_df <- data[, unique(c(id, seq, junc_len, text_fields, num_fields, 
                              seq_fields, cell, locus, traits))]
    tmp_df[[seq]] <- alakazam::maskSeqGaps(tmp_df[[seq]], mask_char=mask_char, 
                                           outer_only=FALSE)
    hc <- dplyr::filter(tmp_df,!!rlang::sym(locus)==rlang::sym(heavy))
    alt <- dplyr::filter(tmp_df,!!rlang::sym(locus)!=rlang::sym(heavy))
    if(nrow(alt) == 0){
      stop(paste("clone",unique(dplyr::pull(data,clone)),
                 "light chain locus not found in dataset!"))
    }
    if(nrow(hc) != 0 & nrow(alt) != 0){
      chain <- "HL"
      print("chain presented as HL. Clone moving forward as if having an HL chain")
    }
    if(length(unique(dplyr::pull(alt,!!locus))) > 1){
        stop(paste("clone",paste(unique(dplyr::pull(data,clone)),collapse=""),
                   "currently only one alternative loci per clone supported"))
      }
  }else{
    # Replace gaps with Ns and masked ragged ends
    tmp_df <- data[, unique(c(id, seq, text_fields, num_fields, seq_fields, traits))]
    tmp_df[[seq]] <- alakazam::maskSeqGaps(tmp_df[[seq]], mask_char=mask_char, 
                                           outer_only=FALSE)
  }
  
  if(chain=="HL"){
    hc[[seq]] <- alakazam::maskSeqEnds(hc[[seq]], mask_char=mask_char,
                                       max_mask=max_mask, trim=FALSE)
    alt[[seq]] <- alakazam::maskSeqEnds(alt[[seq]], mask_char=mask_char, 
                                        max_mask=max_mask, trim=FALSE)
    # Pad ends
    if(pad_end) {
      hc[[seq]] <- alakazam::padSeqEnds(hc[[seq]], pad_char=mask_char, mod3=mod3)
      alt[[seq]] <- alakazam::padSeqEnds(alt[[seq]], pad_char=mask_char, mod3=mod3)
    }
    hc_length <- unique(nchar(dplyr::pull(hc,rlang::sym(seq))))
    alt_length <- unique(nchar(dplyr::pull(alt,rlang::sym(seq))))
    if(length(hc_length) > 1){
      stop(paste("clone",unique(dplyr::pull(data,clone)),
                 "Heavy chain sequences must be same length!"))
    }
    if(length(alt_length) > 1){
      stop(paste("clone",unique(dplyr::pull(data,clone)),
                 "Light chain sequences must be same length!"))
    }
    hc$lsequence <- ""
    hc$hlsequence <- ""
    for(cell_name in unique(dplyr::pull(hc,!!rlang::sym(cell)))){
      if(!cell_name %in% dplyr::pull(alt,rlang::sym(cell))){
        altseq <- paste(rep(mask_char,alt_length),collapse="")
      }else{
        altseq <- dplyr::pull(dplyr::filter(alt,
                                            !!rlang::sym(cell) == cell_name),rlang::sym(seq))
      }
      hc[dplyr::pull(hc,!!rlang::sym(cell)) == cell_name,]$lsequence <- altseq
      hc[dplyr::pull(hc,!!rlang::sym(cell)) == cell_name,]$hlsequence <- 
        paste0(hc[dplyr::pull(hc,!!rlang::sym(cell)) == cell_name,seq],altseq)
    }
    hcd <- dplyr::filter(data,!!rlang::sym(locus)==rlang::sym(heavy))
    altd <- dplyr::filter(data,!!rlang::sym(locus)!=rlang::sym(heavy))
    germline <- alakazam::maskSeqGaps(hcd[[germ]][1], mask_char=mask_char, 
                                      outer_only=FALSE)
    lgermline <- alakazam::maskSeqGaps(altd[[germ]][1], mask_char=mask_char, 
                                       outer_only=FALSE)
    if(pad_end){
      germline <- alakazam::padSeqEnds(germline, pad_char=mask_char, mod3=mod3)
      lgermline <- alakazam::padSeqEnds(lgermline, pad_char=mask_char, mod3=mod3)
    }
    hlgermline <- paste0(germline,lgermline)
    tmp_df <- hc
    loci <- unique(dplyr::pull(data,!!locus))
    tmp_df[[locus]] <- NULL
    tmp_df[[locus]] <- paste(loci,collapse=",")
    chains <- c(rep(unique(dplyr::pull(hc,!!locus)),times=hc_length),
                rep(unique(dplyr::pull(alt,!!locus)),times=alt_length))
    numbers <- c(1:hc_length,1:alt_length)
    if(use_regions){
      hregions <- as.character(
        shazam::setRegionBoundaries(unique(hc[[junc_len]]),
                                    germline,
                                    shazam::IMGT_VDJ_BY_REGIONS)@boundaries)
      lregions <- as.character(
        shazam::setRegionBoundaries(unique(alt[[junc_len]]),
                                    lgermline,
                                    shazam::IMGT_VDJ_BY_REGIONS)@boundaries)
      regions <- c(hregions, lregions)
    } 
    }else if(chain=="L"){
      alt[[seq]] <- alakazam::maskSeqEnds(alt[[seq]], mask_char=mask_char, 
                                          max_mask=max_mask, trim=FALSE)
      # Pad ends
      if(pad_end) {
        #hc[[seq]] <- alakazam::padSeqEnds(hc[[seq]], pad_char=mask_char, mod3=mod3)
        alt[[seq]] <- alakazam::padSeqEnds(alt[[seq]], pad_char=mask_char, mod3=mod3)
      }
      #hc_length <- unique(nchar(dplyr::pull(hc,rlang::sym(seq))))
      alt_length <- unique(nchar(dplyr::pull(alt,rlang::sym(seq))))
      if(length(alt_length) > 1){
        stop(paste("clone",unique(dplyr::pull(data,clone)),
                   "Light chain sequences must be same length!"))
      }
      alt$lsequence <- ""
      alt$hlsequence <- ""
      for(cell_name in unique(dplyr::pull(alt,!!rlang::sym(cell)))){
        if(!cell_name %in% dplyr::pull(alt,rlang::sym(cell))){
          altseq <- paste(rep(mask_char,alt_length),collapse="")
        }else{
          altseq <- dplyr::pull(dplyr::filter(alt,
                                              !!rlang::sym(cell) == cell_name),rlang::sym(seq))
        }
        alt[dplyr::pull(alt,!!rlang::sym(cell)) == cell_name,]$lsequence <- altseq
        alt[dplyr::pull(alt,!!rlang::sym(cell)) == cell_name,]$hlsequence <- altseq
      }
      #hcd <- dplyr::filter(data,!!rlang::sym(locus)==rlang::sym("IGH"))
      altd <- dplyr::filter(data,!!rlang::sym(locus)!=rlang::sym(heavy))
      germline <- ""
      lgermline <- alakazam::maskSeqGaps(altd[[germ]][1], mask_char=mask_char, 
                                         outer_only=FALSE)
      if(pad_end){
        germline <- alakazam::padSeqEnds(germline, pad_char=mask_char, mod3=mod3)
        lgermline <- alakazam::padSeqEnds(lgermline, pad_char=mask_char, mod3=mod3)
      }
      hlgermline <- paste0(germline,lgermline)
      tmp_df <- alt
      loci <- unique(dplyr::pull(data,!!locus))
      tmp_df[[locus]] <- NULL
      tmp_df[[locus]] <- paste(loci,collapse=",")
      chains <- c(rep(unique(dplyr::pull(alt,!!locus)),times=alt_length))
      numbers <- c(1:alt_length)
      if(use_regions){
        hregions <- c()
        lregions <- as.character(
          shazam::setRegionBoundaries(unique(alt[[junc_len]]),
                                      lgermline,
                                      shazam::IMGT_VDJ_BY_REGIONS)@boundaries)
        regions <- c(hregions, lregions)
      }else{
      regions <- rep("N", times=nchar(hlgermline))
      }
      if(length(regions) != nchar(hlgermline)){
      warning(paste("Excluding clone",unique(dplyr::pull(data,clone)),
                    "due to incomplete region definition."))
      return(NULL)
    }
    if(length(chains) != unique(nchar(tmp_df$hlsequence))){
      stop(paste("clone",unique(dplyr::pull(data,clone)),
                 "chains vector not equal to total sequence length!"))
    }
    if(length(chains) != nchar(hlgermline)){
      stop(paste("clone",unique(dplyr::pull(data,clone)),
                 "chains vector not equal to germline sequence length!"))
    }
    new_seq <- "hlsequence"
  }else{
    tmp_df[[seq]] <- alakazam::maskSeqEnds(tmp_df[[seq]], 
                                           mask_char=mask_char, max_mask=max_mask, trim=FALSE)
    if(pad_end){
      tmp_df[[seq]] <- alakazam::padSeqEnds(tmp_df[[seq]], 
                                            pad_char=mask_char, mod3=mod3)
    }
    germline <- alakazam::maskSeqGaps(data[[germ]][1], 
                                      mask_char=mask_char, outer_only=FALSE)
    if(pad_end){
      germline <- alakazam::padSeqEnds(germline, 
                                       pad_char=mask_char, mod3=mod3)
    }
    check <- alakazam::checkColumns(data, c(locus))
    if(check == TRUE){
      loci <- unique(dplyr::pull(data,locus))
      if(length(loci) > 1){
        warning(paste("clone",unique(dplyr::pull(data,clone)),
                      "mutliple loci present but not dealt with!"))
        loci <- paste(loci,collapse=",")
      }
    }else{
      loci <- "N"
    }
    chains <- rep(loci,times=nchar(germline))
    numbers <- 1:nchar(germline) #assumes IMGT numbers
    lgermline <- ""
    hlgermline <- germline
    tmp_df$lsequence <- ""
    tmp_df$hlsequence <- tmp_df[[seq]]
    new_seq <- seq
    if(use_regions){
      regions <- as.character(
        shazam::setRegionBoundaries(unique(data[[junc_len]]),
                                    germline,
                                    shazam::IMGT_VDJ_BY_REGIONS)@boundaries)
    }else{
      regions <- rep("N", times=nchar(germline))
    }
  }
  seq_len <- nchar(tmp_df[[seq]])
  if(any(seq_len != seq_len[1])){
    len_message <- paste0("All sequences are not the same length for data with first ", 
                          id, " = ", tmp_df[[id]][1], ".")
    if (!pad_end){
      len_message <- paste(len_message, 
                           "Consider specifying pad_end=TRUE and verify the multiple alignment.")
    }else{
      len_message <- paste(len_message,
                           "Verify that all sequences are properly multiple-aligned.")
    }
    stop(len_message)
  }
  
  # Remove duplicates
  if(collapse){
    if(is.null(traits)){
      tmp_df <- alakazam::collapseDuplicates(tmp_df, id=id, seq=new_seq, 
                                             text_fields=text_fields, 
                                             num_fields=num_fields, seq_fields=seq_fields,
                                             add_count=add_count, verbose=verbose)
    }else{
      tmp_df <- tmp_df %>%
        dplyr::group_by_at(dplyr::vars(tidyselect::all_of(traits))) %>%
        dplyr::do(alakazam::collapseDuplicates(!!rlang::sym("."), id=id, seq=new_seq, 
                                               text_fields=text_fields, 
                                               num_fields=num_fields, seq_fields=seq_fields,
                                               add_count=add_count, verbose=verbose)) %>%
        dplyr::ungroup()
    }
  }
  if(randomize){
    tmp_df <- tmp_df[sample(1:nrow(tmp_df),replace=FALSE),]
  }
  
  # Define return object
  tmp_names <- names(tmp_df)
  if ("sequence" %in% tmp_names & seq != "sequence") {
    tmp_df <- tmp_df[, tmp_names != "sequence"]
    tmp_names <- names(tmp_df)
  }
  names(tmp_df)[tmp_names == seq] <- "sequence"
  names(tmp_df)[tmp_names == id] <- "sequence_id"
  
  if(chain=="HL"){
    phylo_seq <- "hlsequence"
  }else if(chain=="L"){
    phylo_seq <- "lsequence"
  }else{
    phylo_seq <- "sequence"
  }
  
  if(nrow(tmp_df) == 1 && dup_singles){
    tmp_df2 <- tmp_df
    tmp_df2[[id]] <- paste0(tmp_df[[id]],"_DUPLICATE")
    tmp_df <- bind_rows(tmp_df, tmp_df2)
  }
  
  outclone <- new("airrClone", 
                  data=as.data.frame(tmp_df),
                  clone=as.character(unique(data[[clone]])),
                  germline=alakazam::maskSeqGaps(germline, mask_char=mask_char, 
                                                 outer_only=FALSE),
                  lgermline=alakazam::maskSeqGaps(lgermline, mask_char=mask_char, 
                                                  outer_only=FALSE),
                  hlgermline=alakazam::maskSeqGaps(hlgermline, mask_char=mask_char, 
                                                   outer_only=FALSE), 
                  v_gene=alakazam::getGene(data[[v_call]][1]), 
                  j_gene=alakazam::getGene(data[[j_call]][1]), 
                  junc_len=data[[junc_len]][1],
                  locus=chains,
                  region=regions,
                  numbers=numbers,
                  phylo_seq=phylo_seq)
  
  
  outclone
}

# Remove uniformative columns from data and germline
# 
# \code{cleanAlignment} clean multiple sequence alignments
# @param    clone   \code{airrClone} object
#
# @return   \code{airrClone} object with cleaned alignment
#
cleanAlignment <- function(clone){
  seq <- clone@phylo_seq
  if(seq == "hlsequence"){
    g <- strsplit(clone@hlgermline[1],split="")[[1]]
  }else if(seq == "sequence"){
    g <- strsplit(clone@germline[1],split="")[[1]]
  }else if(seq == "lsequence"){
    g <- strsplit(clone@lgermline[1],split="")[[1]]
  }else{
    stop(paste(seq, "not a recognized sequence type"))
  }
  sk <- strsplit(clone@data[[seq]],split="")
  sites <- seq(1,length(g)-3,by=3)
  ns <- c()
  for(i in sites){ #for each codon site, tally number of NNN codons
    l <- unlist(lapply(sk,function(x) paste(x[i:(i+2)],collapse="")=="NNN"))
    ns <- c(ns,rep(sum(l),length=3))
  }
  # remove uninformative sites from germline and data
  informative <- ns != length(sk)
  l <- lapply(sk,function(x) x=paste(x[informative],collapse=""))
  gm <- paste(g[informative],collapse="")
  
  if(.hasSlot(clone,"locus")){
    clone@locus <- clone@locus[informative]
  }
  if(.hasSlot(clone,"region")){
    clone@region <- clone@region[informative]  
  }
  if(.hasSlot(clone,"numbers")){
    clone@numbers <- clone@numbers[informative]
  }
  if(seq == "hlsequence"){
    clone@hlgermline=gm
  }else if(seq == "sequence"){
    clone@germline=gm
  }else if(seq == "lsequence"){
    clone@lgermline=gm
  }
  clone@data[[seq]]=unlist(l)
  return(clone)
}

#### Preprocessing functions ####

#' Generate an ordered list of airrClone objects for lineage construction
#' 
#' \code{formatClones} takes a \code{data.frame} or \code{tibble} with AIRR or 
#' Change-O style columns as input and masks gap positions, masks ragged ends, 
#' removes duplicates sequences, and merges annotations associated with duplicate
#' sequences. If specified, it will un-merge duplicate sequences with different 
#' values specified in the \code{trait} option. It returns a list of \code{airrClone}
#' objects ordered by number of sequences which serve as input for lineage reconstruction.
#' 
#' @param    data         data.frame containing the AIRR or Change-O data for a clone.
#'                        See \link{makeAirrClone} for required columns and their defaults
#' @param    split_light  split or lump subclones? See \code{getSubclones}.
#' @param    minseq       minimum numbner of sequences per clone
#' @param    majoronly    only return largest subclone and sequences without light chains
#' @param    clone        name of the column containing the identifier for the clone. All 
#'                        entries in this column should be identical.
#' @param    seq          sequence alignment column name.
#' @param    subclone     name of the column containing the identifier for the subclone.
#' @param    chain        if HL, include light chain information if available. If L, only uses 
#'                        light chain information.
#' @param    heavy        name of heavy chain locus (default = "IGH")
#' @param    cell         name of the column containing cell assignment information
#' @param    locus        name of the column containing locus information
#' @param    nproc        number of cores to parallelize formating over.
#' @param    columns      additional data columns to include in output
#' @param    ...          additional arguments to pass to makeAirrClone                        
#'
#' @return   A tibble of \link{airrClone} objects containing modified clones.
#'
#' @details
#' This function is a wrapper for \link{makeAirrClone}. Also removes whitespace,
#' ;, :, and = from ids
#' @seealso  Executes in order \link{makeAirrClone}. Returns a tibble of 
#' \link{airrClone} objects 
#'      which serve as input to \link{getTrees} and \link{findSwitches}.
#' 
#' @examples
#' data(ExampleAirr)
#' # Select two clones, for demonstration purpose
#' sel <- c("3170", "3184")
#' clones <- formatClones(ExampleAirr[ExampleAirr$clone_id %in% sel,],trait="sample_id")
#' @export
formatClones <- function(data, seq="sequence_alignment", clone="clone_id", 
                         subclone="subclone_id",
                         nproc=1, chain="H", heavy="IGH", cell="cell_id", 
                         locus="locus", minseq=2, split_light=FALSE, majoronly=FALSE,
                         columns=NULL, ...) {
  
  if(majoronly){
    if(!subclone %in% names(data)){
      stop("Need subclone designation if majoronly=TRUE")
    }
    data <- filter(data, !!rlang::sym(subclone) <= 1)
  }
  if(chain == "H"){ #if chain is heavy and, discard all non-IGH sequences
    if(!is.null(heavy)){
      if(locus %in% names(data)){
        data <- filter(data, !!rlang::sym(locus) == rlang::sym(heavy))
      }
    }
  }
  if(chain == "HL"){
    if(!subclone %in% names(data)){
      stop("Need subclone designation for heavy+light chain clones")
    }
    if(!locus %in% names(data)){
      stop("Need locus designation for heavy+light chain clones")
    }
    if(is.null(heavy)){
      stop("Need heavy chain (heavy) designation for heavy+light chain clones")
    }    
    lcells <- filter(data,!!rlang::sym(locus)!=rlang::sym(heavy))[[cell]]
    hcells <- filter(data,!!rlang::sym(locus)==rlang::sym(heavy))[[cell]]
    nohcells <- lcells[!lcells %in% hcells]
    if(length(nohcells) > 0){
      data <- filter(data,!(!!rlang::sym(cell) %in% nohcells))
      warning(paste("Removed",length(nohcells),
                    "cells with no heavy chain information"))
    }
  }
  if(chain == "L"){ #if chain is light and, discard all IGH sequences
    if(!is.null(heavy)){
      if(locus %in% names(data)){
        data <- filter(data, !!rlang::sym(locus) != rlang::sym(heavy))
      }
    }
  }
  #edit based on subclone options
  if(split_light){
    if(!subclone %in% names(data)){
      stop("Need subclone designation for heavy+light chain clones")
    }
    if(sum(data[[subclone]] == 0) > 0){
      warning("Assigning subclone 0 (missing light chain) to subclone 1")
      data[data[[subclone]] == 0,][[subclone]] <- 1
    }
    data[[clone]] <- paste0(data[[clone]],"_",data[[subclone]])
  }else if(!split_light && chain=="HL"){
    data <- filter(data, !(!!rlang::sym(locus) != rlang::sym(heavy) &
                             !!rlang::sym(subclone) > 1))
  }
  if(!is.null(columns)){
    if(sum(!columns %in% names(data)) != 0){
      stop(paste("column",
                 paste(columns[!columns %in% names(data)]),
                 "not in data table!"))
    }
  }
  if(sum(is.na(data[[seq]])) > 0){
    warning(paste("Removing",sum(is.na(data[[seq]]))
                  ,"with missing sequences"))
    data <- data[!is.na(data[[seq]]),]
  }
  
  # TODO: Adjust for heavy/light sequences
  counts <- table(data[[clone]])
  rmclones <- names(counts[counts < minseq])
  data <- data[!data[[clone]] %in% rmclones,]
  
  clones <- data %>%
    dplyr::group_by(!!rlang::sym(clone)) %>%
    dplyr::do(data=makeAirrClone(.data, seq=seq,
                                 clone=clone, chain=chain, heavy=heavy, cell=cell, ...))
  saveRDS(clones, file = "~/Downloads/clones.rds")
  
  # remove NULL clone objects
  exclude_clones <- unlist(lapply(clones$data,function(x)is.null(x)))
  if(sum(exclude_clones) > 0){
    warning(paste("Excluding",sum(exclude_clones),"clones"))
  }
  clones <- clones[exclude_clones==F,]
  
  if(nrow(clones) == 0){
    stop("No clones remain after makeAirrClone")
  }
  
  if(chain == "HL"){
    seq_name <- "hlsequence"
  }else if(chain == "H"){
    seq_name <- "sequence"
  }else if(chain == "L"){
    seq_name <- "lsequence"
  }else{
    stop(paste("Chain option",chain,"not recognized"))
  }
  
  fclones <- processClones(clones, nproc=nproc, seq=seq_name, minseq=minseq)
  
  # clone_id must be used for clone ids
  if(clone != "clone_id"){
    fclones$clone_id <- fclones[[clone]]
    fclones <- dplyr::select(fclones, -!!clone)
  }
  
  colpaste <- function(x){
    s <- sort(unique(x))
    if(length(s) > 1){
      paste(s,collapse=",")
    }else{
      s
    }
  }
  if(!is.null(columns)){
    d <- data %>%
      dplyr::select(!!rlang::sym(clone),dplyr::all_of(columns)) %>%
      dplyr::group_by(!!rlang::sym(clone)) %>%
      dplyr::summarize(dplyr::across(dplyr::all_of(columns), dplyr::n_distinct)) %>%
      dplyr::ungroup() %>%
      dplyr::summarize(dplyr::across(dplyr::all_of(columns),max)) %>%
      unlist()
    multi <- names(d[d > 1])
    if(length(multi) > 0){
      warning(paste("columns",paste(multi,collapse=" "),
                    "contain multiple values per clone, flattening with comma"))
    }
    d <- data %>%
      dplyr::select(!!rlang::sym(clone),columns) %>%
      dplyr::group_by(!!rlang::sym(clone)) %>%
      dplyr::summarize(dplyr::across(columns, colpaste))
    
    m <- match(fclones[[clone]],d[[clone]])
    fclones[,columns] <- d[m,columns]
  }
  
  fclones
}


# Clean sequence IDs, order clones by number of sequences, and 
# remove uninformative sites.
# 
# \code{processClones} clean clonal alignments.
# 
# @param    clones   tibble of \code{airrClone} objects containing sequences
# @param    nproc    number of cores for parallelization
# @param    minseq   minimum number of sequences per clone
# @param    seq      column name containing sequence information
#
# @return   a tibble containing \code{airrClone} objects with cleaned sequence
#  aligments
#  
processClones <- function(clones, nproc=1 ,minseq=2, seq){
  
  if(!"tbl" %in% class(clones)){
    print(paste("clones is of class",class(clones)))
    stop("clones must be a tibble of airrClone objects!")
  }else{
    if(class(clones$data[[1]]) != "airrClone"){
      print(paste("clones is list of class",class(clones$data[[1]])))
      stop("clones$data must be a list of airrClone objects!")
    }
  }
  
  threshold <- unlist(lapply(clones$data,function(x)
    length(x@data[[seq]]) >= minseq))
  clones <- clones[threshold,]
  if(nrow(clones) == 0){
    warning(paste("All clones have less than minseq =",minseq,"sequences"))
    return(clones)
  }
  
  clones$data <- lapply(clones$data,function(x){
    x@data$sequence_id=    gsub(":","_",x@data$sequence_id);
    x })
  clones$data <- lapply(clones$data,function(x){
    x@data$sequence_id=gsub(";","_",x@data$sequence_id);
    x })
  clones$data <- lapply(clones$data,function(x){
    x@data$sequence_id=gsub(",","_",x@data$sequence_id);
    x })
  clones$data <- lapply(clones$data,function(x){
    x@data$sequence_id=gsub("=","_",x@data$sequence_id);
    x })
  clones$data <- lapply(clones$data,function(x){
    x@data$sequence_id=gsub(" ","_",x@data$sequence_id);
    x })
  
  max <- max(unlist(lapply(clones$data,function(x)max(nchar(x@data$sequence_id)))))
  if(max > 1000){
    wc <- which.max(unlist(lapply(clones$data,function(x)
      max(nchar(x@data$sequence_id)))))
    stop(paste("Sequence ID of clone",clones$data[[wc]]@clone,"index",
               wc,"too long - over 1000 characters!"))
  }
  
  or <- order(unlist(lapply(clones$data,function(x)nrow(x@data))),
              decreasing=TRUE)
  clones <- clones[or,]
  
  clones$data <- parallel::mclapply(clones$data,
                                    function(x)cleanAlignment(x),mc.cores=nproc)
  
  if(.hasSlot(clones$data[[1]],"locus")){
    clones$locus <- unlist(lapply(clones$data,function(x)
      paste(sort(unique(x@locus)),collapse=",")))
  }
  clones$seqs <- unlist(lapply(clones$data,function(x)nrow(x@data)))
  clones <- dplyr::rowwise(clones)
  clones <- dplyr::ungroup(clones)
  clones
}

# this is to test it 
#fl = formatClones(comb_germline, majoronly=TRUE,verbose=FALSE,
 #                 heavy="IGH", chain="L", cell="cell_id", minseq=1,
  #                columns = c("sample"), dup_singles=TRUE, nproc = 5, 
   #               germ = "germline_alignment_d_mask")

