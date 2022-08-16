# This is the test script to update the makeAirrClone function in dowser. 
# the data used is from the flu data that Ken has access to so I can ensure that we have light chains

library(alakazam)
library(dplyr)
# read in a change-o input and make it light chain 
heavy <- readChangeoDb("~/Documents/Flu Data/heavy.tsv") 
light <- readChangeoDb("~/Documents/Flu Data/light.tsv") 
# make sure there aren't double heavy chains
heavy_count = table(heavy$cell_id)
sum(heavy_count > 1)
multi_heavy_cells = names(heavy_count[heavy_count > 1])
heavy = filter(heavy, !cell_id %in% multi_heavy_cells)

# get the subclones
comb = dowser::getSubclones(heavy,light,cell="cell_id", nproc= 5)
comb = comb[which(!is.na(comb$c_call)),]
# since I already made this bad boy I'm going to down sample this to make it go quicker
to_keep <- sample(comb$cell_id, 500)
comb_sub <- filter(comb, cell_id %in% to_keep)

# get the germlines
references = dowser::readIMGT(dir = "~/share/germlines/imgt/human/vdj")
h = dowser::createGermlines(filter(comb_sub,locus=="IGH"),references, clone = "clone_id", nproc = 5)
k = dowser::createGermlines(filter(comb_sub,locus=="IGK"),references,locus="IGK", clone = "clone_id", nproc = 5)
l = dowser::createGermlines(filter(comb_sub,locus=="IGL"),references,locus="IGL", clone = "clone_id", nproc = 5)
comb_germline = bind_rows(h, l, k)
comb_germline = comb_germline[which(!is.na(comb_germline$c_call)),]


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
data <- filter(data, locus != "IGH")

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
    hc <- dplyr::filter(tmp_df,!!rlang::sym(locus)==rlang::sym("IGH"))
    alt <- dplyr::filter(tmp_df,!!rlang::sym(locus)!=rlang::sym("IGH"))
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
      altd <- dplyr::filter(data,!!rlang::sym(locus)!=rlang::sym("IGH"))
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

test <- makeAirrClone(data, chain = "L")
