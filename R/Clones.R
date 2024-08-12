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
#' @param    v_call       name of the column containing V-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    j_call       name of the column containing J-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    junc_len     name of the column containing the length of the junction as a 
#'                        numeric value. All entries in this column should be identical 
#'                        for any given clone.
#' @param    clone        name of the column containing the identifier for the clone. All 
#'                        entries in this column should be identical.
#' @param    subgroup     name of the column containing the identifier for the subgroup.
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
#' @param    chain        if HL, include light chain information if available.
#' @param    heavy        name of heavy chain locus (default = "IGH")
#' @param    cell         name of the column containing cell assignment information
#' @param    locus        name of the column containing locus information
#' @param    mod3         pad sequences to length mutliple three?
#' @param    randomize    randomize sequence order? Important if using PHYLIP
#' @param    use_regions  assign CDR/FWR regions?
#' @param    dup_singles  Duplicate sequences in singleton clones to include them as trees?
#' @param    light_traits Include the traits from the light chain when concatenating and collapsing trees? 
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
           junc_len="junction_length", clone="clone_id", subgroup = "clone_subgroup",
           mask_char="N", max_mask=0, pad_end=TRUE, text_fields=NULL, num_fields=NULL, 
           seq_fields=NULL, add_count=TRUE, verbose=FALSE, collapse=TRUE, chain="H", 
           heavy=NULL, cell="cell_id", locus="locus", traits=NULL, mod3=TRUE,
           randomize=TRUE, use_regions=TRUE, dup_singles=FALSE, light_traits=FALSE){

    # Check for valid fields
    check <- alakazam::checkColumns(data, 
                                    unique(c(id, seq, germ, v_call, j_call, junc_len, clone, 
                                             text_fields, num_fields, seq_fields, traits)))
    if (check != TRUE) { stop(check) }
    
    if(chain=="HL"){
      check <- alakazam::checkColumns(data, c(cell,locus))
      if (!check) { stop(check) }
      if(is.null(heavy)){
        stop(paste("clone",unique(dplyr::pull(data,clone)),
                   "heavy chain loci ID must be specified if combining loci!"))
      }
      
      heavycount = max(table(data[data[[locus]] == heavy,][[cell]]))
      if(max(heavycount) > 1){
        stop(paste0(sum(heavycount > 1),
                    " cells with multiple heavy chains found. Remove before proceeeding"))
      }
      
      # update the traits and columns CGJ 10/18/23
      if(!is.null(traits) && light_traits){
        new_traits <- c()
        for(i in 1:length(traits)){
          c_trait <- traits[[i]]
          new_traits <- append(new_traits, c_trait)
          new_traits <- append(new_traits, paste0(c_trait, "_light"))
        }
        
        heavy_clones <- dplyr::filter(data,!!rlang::sym(locus)==rlang::sym(heavy))
        light_clones <- dplyr::filter(data,!!rlang::sym(locus)!=rlang::sym(heavy))
        
        for(i in 1:length(traits)){
          heavy_clones[[paste0(traits[i], "_heavy")]] <- heavy_clones[[traits[i]]]
          heavy_clones[[paste0(traits[i], "_light")]] <- unlist(lapply(1:nrow(heavy_clones), function(x){
            cell_id <- heavy_clones[[cell]][x]
            if(cell_id %in% light_clones[[cell]]){
              matching_light <- light_clones[light_clones[[cell]] == cell_id,]
              if(nrow(matching_light) == 1){
                value <- matching_light[[traits[i]]][1]
              } else {
                matching_light <- filter(matching_light, subgroup == min(matching_light[[subgroup]]))
                value <- matching_light[[traits[i]]][1]
              }
            } else {
              value <- NA
            }
            return(value)
          }))
          light_clones[[paste0(traits[i], "_heavy")]] <- NA
          light_clones[[paste0(traits[i], "_light")]] <- light_clones[[traits[i]]]
        }
        
        traits <- new_traits
        data <- rbind(heavy_clones, light_clones)
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
    } else if(chain=="L"){
      check <- alakazam::checkColumns(data, locus)
      if (check != TRUE) { stop(check) }
      
      # Ensure cell and loci columns are not duplicated
      text_fields <- text_fields[text_fields != rlang::sym(cell)]
      text_fields <- text_fields[text_fields != rlang::sym(locus)]
      seq_fields <- seq_fields[seq_fields != rlang::sym(cell)]
      seq_fields <- seq_fields[seq_fields != rlang::sym(locus)]
      # Replace gaps with Ns and masked ragged ends
      tmp_df <- data[, unique(c(id, seq, junc_len, text_fields, num_fields, 
                                seq_fields, locus, traits))]
      tmp_df[[seq]] <- alakazam::maskSeqGaps(tmp_df[[seq]], mask_char=mask_char, 
                                             outer_only=FALSE)
      hc <- dplyr::filter(tmp_df,!!rlang::sym(locus)==rlang::sym(heavy))
      alt <- dplyr::filter(tmp_df,!!rlang::sym(locus)!=rlang::sym(heavy))
      if(nrow(alt) == 0){
        stop(paste("clone",unique(dplyr::pull(data,clone)),
                   "light chain locus not found in dataset!"))
      }
    } else{ 
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
      for(cell_name in 1:nrow(hc)){
        if(is.na(hc[[cell]][cell_name])){
          hc[[cell]][cell_name] <- "bulk"
        }
      }
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
      for(cell_name in 1:nrow(hc)){
        if(hc[[cell]][cell_name] == "bulk"){
          hc[[cell]][cell_name] <- NA
        }
      }
      hcd <- dplyr::filter(data,!!rlang::sym(locus)==rlang::sym(heavy))
      altd <- dplyr::filter(data,!!rlang::sym(locus)!=rlang::sym(heavy))
      
      if(any(hcd[[germ]][1] != hcd[[germ]]) || 
         any(altd[[germ]][1] != altd[[germ]])){
        stop(paste0("Germline sequences for clone ",
                    unique(dplyr::pull(data,clone)),
                    " are not identical. All predicted germline sequences ",
                    "must be identical for each locus within a clone. Be sure to use the ",
                    "createGermlines function before formatClones or makeAirrClone."))
      }
      germline <- alakazam::maskSeqGaps(hcd[[germ]][1], mask_char=mask_char, 
                                        outer_only=FALSE)
      lgermline <- alakazam::maskSeqGaps(altd[[germ]][1], mask_char=mask_char, 
                                         outer_only=FALSE)
      if(pad_end){
        germline <- alakazam::padSeqEnds(germline, pad_char=mask_char, mod3=mod3)
        lgermline <- alakazam::padSeqEnds(lgermline, pad_char=mask_char, mod3=mod3)
        length <- max(c(nchar(germline), max(nchar(hcd[[seq]]))))
        llength <- max(c(nchar(lgermline),max(nchar(altd[[seq]]))))
        if(length > nchar(germline)){
          warning(paste0(
            "Padding germline for clone ",unique(dplyr::pull(data,clone)),
            ", may indicate misalignment.",
            " Should not happen if using createGermlines."))
          germline <- alakazam::padSeqEnds(germline, 
                                           pad_char=mask_char, mod3=mod3, len=length)
        }
        if(llength > nchar(lgermline)){
          warning(paste0(
            "Padding germline for clone ",unique(dplyr::pull(data,clone)),
            ", may indicate misalignment.",
            "Should not happen if using createGermlines."))
          lgermline <- alakazam::padSeqEnds(lgermline, 
                                            pad_char=mask_char, mod3=mod3, len=llength)
        }
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
    } else if(chain=="L"){
      tmp_df[[seq]] <- alakazam::maskSeqEnds(tmp_df[[seq]], 
                                             mask_char=mask_char, max_mask=max_mask, trim=FALSE)
      if(pad_end){
        tmp_df[[seq]] <- alakazam::padSeqEnds(tmp_df[[seq]], 
                                              pad_char=mask_char, mod3=mod3)
      }
      if(any(data[[germ]][1] != data[[germ]])){
        stop(paste0("Germline sequences for clone ",
                    unique(dplyr::pull(data,clone)),
                    " are not identical. All predicted germline sequences ",
                    "must be identical within a clone. Be sure to use the ",
                    "createGermlines function before formatClones or makeAirrClone."))
      }
      lgermline <- alakazam::maskSeqGaps(data[[germ]][1], 
                                         mask_char=mask_char, outer_only=FALSE)
      if(pad_end){
        lgermline <- alakazam::padSeqEnds(lgermline, 
                                          pad_char=mask_char, mod3=mod3)
        length <- max(c(nchar(lgermline),max(nchar(tmp_df[[seq]]))))
        if(length > nchar(lgermline)){
          warning(paste0(
            "Padding germline for clone ",unique(dplyr::pull(data,clone)),
            ", may indicate misalignment.",
            " Should not happen if using createGermlines."))
          germline <- alakazam::padSeqEnds(lgermline, 
                                           pad_char=mask_char, mod3=mod3, len=length)
        }
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
      chains <- rep(loci,times=nchar(lgermline))
      numbers <- 1:nchar(lgermline) #assumes IMGT numbers
      germline <- ""
      hlgermline <- lgermline
      tmp_df$lsequence <- tmp_df[[seq]]
      tmp_df$hlsequence <- tmp_df[[seq]]
      new_seq <- seq
      if(use_regions){
        regions <- as.character(
          shazam::setRegionBoundaries(unique(data[[junc_len]]),
                                      lgermline,
                                      shazam::IMGT_VDJ_BY_REGIONS)@boundaries)
      }else{
        regions <- rep("N", times=nchar(lgermline))
      }
    }else{ 
      tmp_df[[seq]] <- alakazam::maskSeqEnds(tmp_df[[seq]], 
                                             mask_char=mask_char, max_mask=max_mask, trim=FALSE)
      if(pad_end){
        tmp_df[[seq]] <- alakazam::padSeqEnds(tmp_df[[seq]], 
                                              pad_char=mask_char, mod3=mod3)
      }
      if(any(data[[germ]][1] != data[[germ]])){
        stop(paste0("Germline sequences for clone ",
                    unique(dplyr::pull(data,clone)),
                    " are not identical. All predicted germline sequences ",
                    "must be identical within a clone. Be sure to use the ",
                    "createGermlines function before formatClones or makeAirrClone."))
      }
      germline <- alakazam::maskSeqGaps(data[[germ]][1], 
                                        mask_char=mask_char, outer_only=FALSE)
      if(pad_end){
        germline <- alakazam::padSeqEnds(germline, 
                                         pad_char=mask_char, mod3=mod3)
        length <- max(c(nchar(germline),max(nchar(tmp_df[[seq]]))))
        if(length > nchar(germline)){
          warning(paste0(
            "Padding germline for clone ",unique(dplyr::pull(data,clone)),
            ", may indicate misalignment.",
            " Should not happen if using createGermlines."))
          germline <- alakazam::padSeqEnds(germline, 
                                           pad_char=mask_char, mod3=mod3, len=length)
        }
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
      tmp_df <- dplyr::bind_rows(tmp_df, tmp_df2)
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
    lg <- strsplit(clone@lgermline[1],split="")[[1]]
    hg <- strsplit(clone@germline[1],split="")[[1]]
  }else if(seq == "sequence"){
    g <- strsplit(clone@germline[1],split="")[[1]]
  }else if(seq == "lsequence"){
    g <- strsplit(clone@lgermline[1],split="")[[1]]
  }else{
    stop(paste(seq, "not a recognized sequence type"))
  }
  if(seq == "hlsequence"){
    sk <- strsplit(clone@data[[seq]],split="")
    sites <- seq(1,length(g)-3,by=3)
    sk_h <- strsplit(clone@data$sequence,split="")
    sk_l <- strsplit(clone@data$lsequence,split="")
    h_sites <- seq(1,length(hg)-3,by=3)
    l_sites <- seq(1,length(lg)-3,by=3)
    ns <- c()
    for(i in sites){ #for each codon site, tally number of NNN codons
      l <- unlist(lapply(sk,function(x) paste(x[i:(i+2)],collapse="")=="NNN"))
      ns <- c(ns,rep(sum(l),length=3))
    }
    hns <- c()
    for(i in h_sites){ #for each codon site, tally number of NNN codons
      l <- unlist(lapply(sk_h,function(x) paste(x[i:(i+2)],collapse="")=="NNN"))
      hns <- c(hns,rep(sum(l),length=3))
    }
    lns <- c()
    for(i in l_sites){ #for each codon site, tally number of NNN codons
      l <- unlist(lapply(sk_l,function(x) paste(x[i:(i+2)],collapse="")=="NNN"))
      lns <- c(lns,rep(sum(l),length=3))
    }
  } else{
    sk <- strsplit(clone@data[[seq]],split="")
    sites <- seq(1,length(g)-3,by=3)
    ns <- c()
    for(i in sites){ #for each codon site, tally number of NNN codons
      l <- unlist(lapply(sk,function(x) paste(x[i:(i+2)],collapse="")=="NNN"))
      ns <- c(ns,rep(sum(l),length=3))
    }
  }

  # remove uninformative sites from germline and data
  informative <- ns != length(sk)
  l <- lapply(sk,function(x) x=paste(x[informative],collapse=""))
  gm <- paste(g[informative],collapse="")
  if(seq == "hlsequence"){
    hinformative <- hns != length(sk_h)
    hl <- lapply(sk_h,function(x) x=paste(x[hinformative],collapse=""))
    hgm <- paste(hg[hinformative],collapse="")
    linformative <- lns != length(sk_l)
    ll <- lapply(sk_l,function(x) x=paste(x[linformative],collapse=""))
    lgm <- paste(lg[linformative],collapse="")
  }
  
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
    clone@lgermline=lgm
    clone@germline=hgm
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
#' values specified in the \code{traits} option. It returns a list of \code{airrClone}
#' objects ordered by number of sequences which serve as input for lineage reconstruction.
#' 
#' @param    data         data.frame containing the AIRR or Change-O data for a clone.
#'                        See \link{makeAirrClone} for required columns and their defaults
#' @param    split_light  split or lump subgroups? See \code{resolveLightChains}.
#' @param    filterstop   only use sequences that do not contain an in-frame stop codon
#' @param    minseq       minimum number of sequences per clone
#' @param    majoronly    only return largest subgroup and sequences without light chains
#' @param    clone        name of the column containing the identifier for the clone. All 
#'                        entries in this column should be identical.
#' @param    seq          sequence alignment column name.
#' @param    subgroup     name of the column containing the identifier for the subgroup.
#' @param    chain        if HL, include light chain information if available.
#' @param    heavy        name of heavy chain locus (default = "IGH")
#' @param    cell         name of the column containing cell assignment information
#' @param    locus        name of the column containing locus information
#' @param    nproc        number of cores to parallelize formating over.
#' @param    columns      additional data columns to include in output
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
#' @param    chain        if HL, include light chain information if available.
#' @param    heavy        name of heavy chain locus (default = "IGH")
#' @param    mod3         pad sequences to length mutliple three?
#' @param    randomize    randomize sequence order? Important if using PHYLIP
#' @param    use_regions   assign CDR/FWR regions?
#' @param    dup_singles   Duplicate sequences in singleton clones to include them as trees?   
#' @param    light_traits  Include the traits from the light chain when concatenating and collapsing trees?              
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
#' clones <- formatClones(ExampleAirr[ExampleAirr$clone_id %in% sel,],traits="sample_id")
#' @export
formatClones <- function(data, seq="sequence_alignment", clone="clone_id", 
    subgroup="clone_subgroup", id="sequence_id", 
    germ="germline_alignment_d_mask", v_call="v_call", j_call="j_call",
    junc_len="junction_length", mask_char="N",
    max_mask=0, pad_end=TRUE, text_fields=NULL, num_fields=NULL, seq_fields=NULL,
    add_count=TRUE, verbose=FALSE, collapse=TRUE,
    cell="cell_id", locus="locus", traits=NULL, mod3=TRUE, randomize=TRUE,
    use_regions=TRUE, dup_singles=FALSE, nproc=1, chain="H", heavy="IGH", 
    filterstop=FALSE, minseq=2, split_light=FALSE, light_traits=FALSE, majoronly=FALSE,
    columns=NULL){
  
  if(majoronly){
    if(!subgroup %in% names(data)){
      stop("Need subgroup designation if majoronly=TRUE")
    }
    data <- dplyr::filter(data, !!rlang::sym(subgroup) <= 1)
  }
  if(filterstop){
    full_nrow <- nrow(data)
    data <- parallel::mclapply(1:nrow(data), function(x){
      sub_seq <- alakazam::translateDNA(data[[seq]][x])
      if(!grepl("\\*", sub_seq)){
        data[x,] 
      }
    }, mc.cores = nproc)
    data <- do.call(rbind, data)
    if(nrow(data) != full_nrow){
      n_removed <- full_nrow - nrow(data)
      warning(paste0("There were ", n_removed, " sequence(s) with an inframe stop codon",
      " and were removed. If you want to keep these sequences use the option filterstop=FALSE."))
    }
  }
  if(!clone %in% names(data)){
    stop(clone," column not found.")
  }
  
  # CGJ 8/10/23
  # added factor check for trait and fields columns 
  if(!is.null(traits) || !is.null(text_fields) || !is.null(num_fields) || !is.null(seq_fields)){
    if(!is.null(traits)){
      sub_data <- data[, traits]
      check <- sapply(sub_data, is.factor)
      if(TRUE %in% check){
        factor_traits <- traits[check]
        stop(paste("Traits cannot be factors. Some indicated trait variable(s):", factor_traits, "have been detected to be factors."))
      }
    }
    if(!is.null(text_fields)){
      sub_data <- data[, text_fields]
      check <- sapply(sub_data, is.factor)
      if(TRUE %in% check){
        factor_text <- text_fields[check]
        stop(paste("text_fields cannot be factors. Some indicated text_field variable(s):", factor_text, "have been detected to be factors."))
      }
    }
    if(!is.null(num_fields)){
      sub_data <- data[, num_fields]
      check <- sapply(sub_data, is.factor)
      if(TRUE %in% check){
        factor_num <- num_fields[check]
        stop(paste("num_fields cannot be factors. Some indicated num_fields variable(s):", factor_num, "have been detected to be factors."))
      }
    }
    if(!is.null(seq_fields)){
      sub_data <- data[, seq_fields]
      check <- sapply(sub_data, is.factor)
      if(TRUE %in% check){
        factor_seq <- seq_fields[check]
        stop(paste("seq_fields cannot be factors. Some indicated seq_fields variable(s):", factor_seq, "have been detected to be factors."))
      }
    }
  }

  # CGJ 7/25/23 changed ordering and moved away from dplyr filtering step 
  # caused a segfault due to 'memory not mapped'
  # base R's filtering doesn't do this. 
  if(chain == "H"){ #if chain is heavy and, discard all non-IGH sequences
    if(!is.null(heavy)){
      if(locus %in% names(data)){
#        data <- dplyr::filter(data, !!rlang::sym(locus) == rlang::sym(heavy))
        data <- data[data[[locus]] == rlang::sym(heavy),]
      }
    }
  }
  if(chain == "HL"){
    if(!subgroup %in% names(data)){
      stop("Need subgroup designation for heavy+light chain clones")
    }
    if(!locus %in% names(data)){
      stop("Need locus designation for heavy+light chain clones")
    }
    if(is.null(heavy)){
      stop("Need heavy chain (heavy) designation for heavy+light chain clones")
    }    
    lcells <- dplyr::filter(data,!!rlang::sym(locus)!=rlang::sym(heavy))[[cell]]
    hcells <- dplyr::filter(data,!!rlang::sym(locus)==rlang::sym(heavy))[[cell]]
    nohcells <- lcells[!lcells %in% hcells]
    if(length(nohcells) > 0){
      data <- dplyr::filter(data,!(!!rlang::sym(cell) %in% nohcells))
      warning(paste("Removed",length(nohcells),
                    "cells with no heavy chain information"))
    }
  }
  if(chain == "L"){ #if chain is light and, discard all IGH sequences
    if(!is.null(heavy)){
      if(locus %in% names(data)){
        data <- dplyr::filter(data, !!rlang::sym(locus) != rlang::sym(heavy))
      }
    }
  }
  #edit based on subgroup options
  if(split_light){
    if(!subgroup %in% names(data)){
      stop("Need subgroup designation for heavy+light chain clones")
    }
    if(sum(data[[subgroup]] == 0) > 0){
      warning("Assigning subgroup 0 (missing light chain) to subgroup 1")
      data[data[[subgroup]] == 0,][[subgroup]] <- 1
    }
    data[[clone]] <- paste0(data[[clone]],"_",data[[subgroup]])
  }else if(!split_light && chain=="HL"){
    #since we're not splitting the light chains, can only include the biggest subgroup
    #for tree building
    data <- dplyr::filter(data, !(!!rlang::sym(locus) != rlang::sym(heavy) &
                                    !!rlang::sym(subgroup) > 1))
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
    dplyr::do(data=makeAirrClone(.data, seq=seq, germ=germ, id=id,
      v_call=v_call, j_call=j_call, junc_len=junc_len, mask_char=mask_char,
      max_mask=max_mask, pad_end=pad_end, text_fields=text_fields, num_fields=num_fields,
      seq_fields=seq_fields, add_count=add_count, verbose=verbose, collapse=collapse,
      heavy=heavy, cell=cell, locus=locus, traits=traits, mod3=mod3, randomize=randomize,
      use_regions=use_regions, dup_singles=dup_singles,
      clone=clone, chain=chain))
  
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
                    "contain multiple values per clone, setting to character, flattening with comma"))
      for(col in columns){
        data[[col]] = as.character(data[[col]])
      }
    }
    d <- data %>%
      dplyr::select(!!rlang::sym(clone),dplyr::all_of(columns)) %>%
      dplyr::group_by(!!rlang::sym(clone)) %>%
      dplyr::summarize(dplyr::across(columns, colpaste))
    
    m <- match(fclones$clone_id,d[[clone]])
    fclones[,columns] <- d[m,columns]
  }
  
  fclones
}


#' \code{maskCodons} Masks codons split by insertions
#' @param    id              sequence id
#' @param    q               (query) un-aligned input sequence (sequence)
#' @param    s               (subject) aligned input sequence (sequence_alignment)
#' @param    keep_alignment  store q and s alignments
#' @param    keep_insertions return removed insertion sequences?
#' @param    gap_opening      gap opening penalty (Biostrings::pairwiseAlignment)
#' @param    gap_extension    gap extension penalty (Biostrings::pairwiseAlignment)
#' @param    mask            if FALSE, don't mask codons
#' @return   A list with split codons masked, if found (sequence_masked).
#'
#' @details
#' Performs global alignment of q and s, masks codons in s that are split by 
#' insertions (see example)
#' masking_note notes codon positions in subject_alignment sequence that were 
#' masked, if found.
#' subject_alignment contains subject sequence aligned to query (q) sequence
#' query_alignment contains query sequence aligned to subject (q) sequence
#' sequence_masked will be NA if frameshift or alignment error detected/
#' @seealso  \link{maskSequences}, Biostrings::pairwiseAlignment.
#' 
#' @examples
#' s = "ATCATCATC..."
#' q = "ATCTTTATCATC"
#' print(maskCodons(1,q,s,TRUE))
#'
#' s <- "ATCATCATC..."
#' q <- "ATTTTCATCATC"
#' print(maskCodons("test",q,s,keep_alignment=TRUE,keep_insertions=TRUE))
#' @export
maskCodons <- function(id, q, s, keep_alignment=FALSE, gap_opening=5, 
                       gap_extension=1, keep_insertions=FALSE, mask=TRUE){
  results <- list(
    sequence_id =id,
    sequence_masked="",
    masking_note="",
    insertions="",
    subject_alignment="",
    query_alignment="")
  
  # remove in-frame IMGT gaps
  sg <- gsub("\\.\\.\\.","",s)
  
  # if sequences are identical, nothing to fix
  if(sg == q || !mask){
    results$subject_alignment <- sg
    results$query_alignment <- q
    results$sequence_masked_v <- s
    return(results)
  }
  
  # store in-frame IMGT gap positions
  gaps <- stringr::str_locate_all(s,"\\.\\.\\.")
  
  # remove in-frame gaps
  sgf <- gsub("---", "", sg)
  if(grepl("-",sgf)){ # if read-shifting gaps present, quit
    results$sequence_masked_v <- NA
    results$masking_note <- "Frameshift in sequence"
    return(results)
  }
  # mark in-frame gaps to keep them during alignment
  sg <- gsub("---", "XXX", sg)
  
  # perform global alignment
  n <- Biostrings::pairwiseAlignment(q, sg, type="global",
                                     gapOpening=gap_opening, gapExtension=gap_extension)
  qa <- as.character(n@pattern)
  sa <- as.character(n@subject)
  if(keep_alignment){
    results$subject_alignment <- sa
    results$query_alignment <- qa
  }
  
  if(keep_insertions){
    insertions <- stringr::str_locate_all(sa,"\\-+")[[1]]
    indels <- ""
    if(nrow(insertions) > 0){
      for(i in 1:nrow(insertions)){
        ins <- substr(qa,insertions[i,1],insertions[i,2])
        if(i == 1){
          indels <- paste0(insertions[i,1],"-",ins)
        }else{
          indels <- paste0(indels,",",
                           paste0(insertions[i,1],"-",ins))
        }
      }
    }
    results$insertions <- indels
  }
  
  # if s began out of frame, add IMGT dots back.
  if(grepl("^\\.\\.",sg)){
    sa <- paste0("..",sa)
    if(keep_alignment){
      results$subject_alignment <- sa
    }
  }else if(grepl("^\\.",sg)){
    sa <- paste0(".",sa)
    if(keep_alignment){
      results$subject_alignment <- sa
    }
  }
  
  # if alignment is poor, pairwiseAlignment will trim sequences
  # freak out and die if this happens
  if(nchar(sa) < nchar(sg) || nchar(qa) < nchar(q)){
    results$sequence_masked_v <- NA
    results$masking_note <- "Alignment error"
    return(results)
  }
  
  # check for frameshifts after alignment
  sgf <- gsub("---", "", sa)
  if(grepl("-",sgf)){
    results$sequence_masked_v <- NA
    results$masking_note <- "Frameshift after alignment"
    return(results)
  }
  
  # if aligned sequences differ by a gap character, maybe mask
  if(qa != sa && grepl("\\-",sa)){
    sas <- strsplit(sa,split="")[[1]]
    mask <- c()
    nseq <- c()
    # loop through subject alignment one codon at a time.
    for(j in 1:ceiling(length(sas)/3)){
      index <- (j-1)*3 + 1
      triple <- paste0(sas[index:(index+2)],collapse="")
      triple <- gsub("NA","",triple) # happens if sequence isn't mod 3
      # if codon contains > 0 but < 3 gap characters, mask it.
      m  <- 0
      if(grepl("\\-",triple)){
        triple <- gsub("[A-Z]","N",triple) #mask characters
        triple <- gsub("\\-","",triple) #discard gaps
        if(nchar(triple) != 0){ #if anything left, it's a masked codon
          m <- 1
        }
      }
      mask <- c(mask,m)
      nseq <- c(nseq,triple) #build new sequence
    }
    maskseq <- paste0(nseq,collapse="")
    
    # masked sequence should be same length as sg
    if(nchar(maskseq) != nchar(sg)){
      print(paste(maskseq,"\n",sg))
      stop("Sequence masking failed")
    }
    
    # add IMGT gaps back in
    sequence_alignment <- maskseq
    if(nrow(gaps[[1]]) > 0){
      for(j in 1:nrow(gaps[[1]])){
        sequence_alignment <- 
          paste0(substr(sequence_alignment,1,gaps[[1]][j,1]-1),
                 "...",substr(sequence_alignment,gaps[[1]][j,1],
                              nchar(sequence_alignment)))
      }
    }
    
    # convert marked gap characters back into gaps
    sequence_alignment <- gsub("X","-",sequence_alignment)
    
    # masked sequence should have no mismatched from starting sequence
    if(alakazam::seqDist(sequence_alignment,s) != 0){
      print(paste(sequence_alignment,"\n",s))
      stop("Adding gaps failed")
    }
    results$sequence_masked_v <- sequence_alignment
    if(sum(mask == 1) > 0){ #note which positions were masked
      results$masking_note <- 
        paste(which(mask == 1),collapse=",")
    }
    return(results)
  }else{
    results$sequence_masked_v <- s
    return(results)
  }
}

#' \code{maskSequences} Mask codons split by insertions in V gene
#' @param data                BCR data table
#' @param sequence_id         sequence id column
#' @param sequence            input sequence column (query)
#' @param sequence_alignment  aligned (IMGT-gapped) sequence column (subject)
#' @param v_sequence_start    V gene start position in sequence
#' @param v_sequence_end      V gene end position in sequence
#' @param v_germline_start    V gene start position in sequence_alignment
#' @param v_germline_end      V gene end position in sequence_alignment
#' @param keep_alignment      store alignment of query and subject sequences?
#' @param keep_insertions     return removed insertion sequences?
#' @param mask_codons         mask split codons?
#' @param mask_cdr3           mask CDR3 sequences?
#' @param junction_length     name of junction_length column
#' @param nproc               number of cores to use
#' @return   A tibble with masked sequence in sequence_masked column, 
#'  as well as other columns.
#'
#' @details
#' Performs global alignment of sequence and sequence_alignment, 
#' masking codons in sequence_alignment that are split by insertions (see examples)
#' masking_note notes codon positions in subject_alignment sequence that 
#' were masked, if found.
#' subject_alignment contains subject sequence aligned to query sequence (only 
#' if keep_alignment=TRUE)
#' query_alignment contains query sequence aligned to subject sequence (only if 
#' keep_alignment=TRUE)
#' sequence_masked will be NA if frameshift or alignment error detected. This 
#' will be noted
#' insertions column will be returned if keep_insertions=TRUE, contains a
#' comma-separated list of each <position in query alignment>-<sequence>. See example.
#' in masking_note.
#' @seealso  \link{maskCodons}, Biostrings::pairwiseAlignment.
#' 
#' @export
maskSequences <- function(data,  sequence_id = "sequence_id", sequence = "sequence",
                          sequence_alignment="sequence_alignment", v_sequence_start = "v_sequence_start", 
                          v_sequence_end = "v_sequence_end", v_germline_start = "v_germline_start", 
                          v_germline_end = "v_germline_end", junction_length="junction_length",
                          keep_alignment = FALSE, keep_insertions=FALSE, 
                          mask_codons=TRUE, mask_cdr3=TRUE, nproc=1){
  
  ids <- data[[sequence_id]]
  qi <- substr(data[[sequence]],
               data[[v_sequence_start]],data[[v_sequence_end]])
  # v segment
  si <- substr(data[[sequence_alignment]],
               data[[v_germline_start]],
               data[[v_germline_end]])
  # rest of the sequence
  ei <- substr(data[[sequence_alignment]],
               data[[v_germline_end]]+1,
               nchar(data[[sequence_alignment]]))
  
  if(max(table(ids)) > 1){
    stop("Sequence IDs are not unique")
  }
  
  results <- dplyr::bind_rows(
    parallel::mclapply(1:length(qi),function(x){
      mask <- maskCodons(ids[x], qi[x], si[x],
                         keep_alignment=keep_alignment,
                         keep_insertions=keep_insertions,
                         mask=mask_codons)
      if(is.na(mask$sequence_masked_v)){
        mask$sequence_masked <- NA
      }else{
        mask$sequence_masked <- 
          paste0(mask$sequence_masked_v,ei[x])
      }
      mask
    }, mc.cores=nproc))
  
  m <- match(data[[sequence_id]], results[[sequence_id]])
  
  if(sum(results[m,]$sequence_id != data$sequence_id) > 0){
    stop("Sequence ids don't match")
  }
  
  data <- bind_cols(data,
                    sequence_masked=results[m,]$sequence_masked,
                    masking_note=results [m,]$masking_note)
  if(keep_alignment){
    data <- bind_cols(data,
                      subject_alignment=results[m,]$subject_alignment,
                      query_alignment=results[m,]$query_alignment)
  }
  if(keep_insertions){
    data <- bind_cols(data,
                      insertions=results[m,]$insertions)
  }
  
  if(mask_cdr3){
    data$sequence_masked <- unlist(lapply(1:nrow(data),function(x){
      if(is.na(data$sequence_masked[x])){
        return(data$sequence_masked[x])
      }
      regions <- as.character(
        shazam::setRegionBoundaries(data[[junction_length]][x],
                                    data$sequence_masked[x],
                                    shazam::IMGT_VDJ_BY_REGIONS)@boundaries)
      if(!is.na(data$sequence_masked[x]) && sum(regions == "cdr3") > 0){
        s <- strsplit(data$sequence_masked[x],split="")[[1]]
        s[regions == "cdr3"] = "N"
        s <- paste(s, collapse="")
      }else{
        s <- data$sequence_masked[x]
      }
      s
    }))
  }
  
  include <- !is.na(data$sequence_masked)
  if(sum(include) == 0){
    warning("Masking failed for all sequences")
    return(data)
  }
  
  #masked sequences should be same length as sequence_alignment
  diffs <- nchar(data$sequence_alignment[include]) - 
    nchar(data$sequence_masked[include])
  if(sum(diffs) > 0){
    print(data[diffs > 0,]$sequence_id)
    stop("Error in masking above sequences (length)")
  }
  
  #masked sequences should have no mismatches from sequence_alignment
  dists <- unlist(parallel::mclapply(1:nrow(data[include,]), function(x)
    alakazam::seqDist(data$sequence_alignment[include][x],
                      data$sequence_masked[include][x]),mc.cores=nproc))
  if(sum(dists) > 0){
    print(data[dists > 0,]$sequence_id)
    stop("Error in masking above sequences (mismatches)")
  }
  return(data)
}

#' #' Deprecated! Use resolveLightChains
#' 
#' \code{getSubClones} plots a tree or group of trees
#' @param    heavy        a tibble containing heavy chain sequences with clone_id
#' @param    light        a tibble containing light chain sequences
#' @param    nproc        number of cores for parallelization
#' @param    minseq       minimum number of sequences per clone
#' @param    id           name of the column containing sequence identifiers.
#' @param    seq          name of the column containing observed DNA sequences. All 
#'                        sequences in this column must be multiple aligned.
#' @param    clone        name of the column containing the identifier for the clone. All 
#'                        entries in this column should be identical.
#' @param    cell         name of the column containing identifier for cells.
#' @param    v_call       name of the column containing V-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    j_call       name of the column containing J-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    junc_len     name of the column containing the length of the junction as a 
#'                        numeric value. All entries in this column should be identical 
#'                        for any given clone.
#' @param    nolight      string to use to indicate a missing light chain
#'
#' @return   a tibble containing 

#' @export

getSubclones <- function(heavy, light, nproc=1, minseq=1,
                         id="sequence_id", seq="sequence_alignment",
                         clone="clone_id", cell="cell_id", v_call="v_call", j_call="j_call",
                         junc_len="junction_length", nolight="missing"){
  stop("This function has been depreciated. Please use resolveLightChains.")
}


#' Define subgroups within clones based on light chain rearrangements
#' 
#' \code{resolveLightChains} resolve light chain V and J subgroups within a clone
#' @param    data         a tibble containing heavy and light chain sequences with clone_id
#' @param    nproc        number of cores for parallelization
#' @param    minseq       minimum number of sequences per clone
#' @param    locus        name of column containing locus values
#' @param    heavy        value of heavy chains in locus column. All other values will be 
#'                        treated as light chains
#' @param    id           name of the column containing sequence identifiers.
#' @param    seq          name of the column containing observed DNA sequences. All 
#'                        sequences in this column must be multiple aligned.
#' @param    clone        name of the column containing the identifier for the clone. All 
#'                        entries in this column should be identical.
#' @param    cell         name of the column containing identifier for cells.
#' @param    v_call       name of the column containing V-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    j_call       name of the column containing J-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    junc_len     name of the column containing the length of the junction as a 
#'                        numeric value. All entries in this column should be identical 
#'                        for any given clone.
#' @param    nolight      string to use to indicate a missing light chain
#' @param    pad_ends          pad sequences within a clone to same length?
#'
#' @return   a tibble containing the same data as inputting, but with the column clone_subgroup
#' added. This column contains subgroups within clones that contain distinct light chain
#' V and J genes, with at most one light chain per cell.
#' @details
#' 1. Make temporary array containing light chain clones
#' 2. Enumerate all possible V, J, and junction length combinations
#' 3. Determine which combination is the most frequent
#' 4. Assign sequences with that combination to clone t
#' 5. Copy those sequences to return array
#' 6. Remove all cells with that combination from temp array
#' 7. Repeat 1-6 until temporary array zero.
#' If there is more than rearrangement with the same V/J
#' in the same cell, pick the one with the highest non-ambiguous
#' characters. Cells with missing light chains are grouped with their
#' subgroup with the closest matching heavy chain (Hamming distance)
#' then the largest and lowest index subgroup if ties are present.
#' 
#' Outputs of the function are 
#' 1. clone_subgroup which identifies the light chain VJ rearrangement that sequence belongs to within it's clone
#' 2. clone_subgroup_id which combines the clone_id variable and the clone_subgroup variable by a "_". 
#' 3. vj_cell which combines the vj_gene and vj_alt_cell columns by a ",".
# TODO: add "fields" option consistent with other functions
#' @export
resolveLightChains <- function(data, nproc=1, minseq=1,locus="locus",heavy="IGH",
                               id="sequence_id", seq="sequence_alignment",
                               clone="clone_id", cell="cell_id", v_call="v_call", j_call="j_call",
                               junc_len="junction_length", nolight="missing", pad_ends=TRUE){
  
  subgroup <- "clone_subgroup"

  light <- data[data[[locus]] != heavy,]
  heavy <- data[data[[locus]] == heavy,]

  if(nrow(heavy) == 0){
    stop("No heavy chains found in data")
  }
  if(nrow(light) == 0){
    warning("No light chains found in data! Assigning all sequences to subgroup 1.")
    data[[subgroup]] = 1
    return(data)
  }

  scount <- table(heavy[[clone]])
  big <- names(scount)[scount >= minseq]
  heavy <- dplyr::filter(heavy,(!!rlang::sym(clone) %in% big))
  
  if(max(table(heavy[[id]])) > 1){
    stop("Sequence IDs in heavy dataframe must be unique!")
  }
  if(max(table(light[[id]])) > 1){
    stop("Sequence IDs in light dataframe must be unique!")
  }
  
  heavycount = table(heavy[[cell]])
  if(max(heavycount) > 1){
    stop(paste0(sum(heavycount > 1),
                " cells with multiple heavy chains found. Remove before proceeeding"))
  }
  
  # filter out light chains with missing cell IDs
  missing_cell <- light[is.na(light[[cell]]),]
  if(nrow(missing_cell) > 0){
    warning(paste("removing",nrow(missing_cell),"light chains with missing cell IDs."))
    light <- light[!is.na(light[[cell]]),]
  }

  heavy$vj_gene <- nolight
  heavy$vj_alt_cell <- NA # set these to NA, since they're pretty rare
  heavy[[subgroup]] <- 1
  light$vj_gene <- nolight
  light$vj_alt_cell <- NA
  light[[subgroup]] <- 1
  light[[clone]] <- -1
  paired <- parallel::mclapply(unique(heavy[[clone]]),function(cloneid){
    # Get heavy chains within a clone, and corresponding light chains
    # separate heavy chains with (sc) and without (bulk) paired light chains
    hd <- dplyr::filter(heavy,!!rlang::sym(clone) == cloneid)
    ld <- dplyr::filter(light,!!rlang::sym(cell) %in% hd[[!!cell]])
    ld <- dplyr::filter(ld, !is.na(!!rlang::sym(cell)))
    hd_sc <- hd[hd[[cell]] %in% ld[[cell]] & !is.na(hd[[cell]]),] # added is.na(cell) catch
    hd_bulk <- hd[!hd[[cell]] %in% ld[[cell]] | is.na(hd[[cell]]),]
    if(nrow(ld) == 0){
      hd$clone_subgroup_id <- paste0(hd[[clone]],"_",hd[[subgroup]])
      hd$vj_cell <- sapply(1:nrow(hd), function(x){
        if(!is.na(hd$vj_alt_cell[x])){
          paste(hd$vj_gene[x],hd$vj_alt_cell[x],sep=",")        
        }else{
          hd$vj_gene[x]
        }
      })
      return(hd)
    }
    ltemp <- dplyr::filter(ld, !is.na(!!rlang::sym(cell)))
    ltemp[[clone]] <- -1
    ld <- dplyr::tibble()
    lclone <- 1
    while(nrow(ltemp) > 0){
      #expand ambiguous V/J calls
      lvs <- strsplit(ltemp[[v_call]],split=",")
      ljs <- strsplit(ltemp[[j_call]],split=",")
      jlens <- ltemp[[junc_len]]
      # get all combinations of V/J calls for each light chain
      combos <-
        lapply(1:length(lvs),function(w)
          unlist(lapply(lvs[[w]],function(x)
            unlist(lapply(ljs[[w]],function(y)
              lapply(jlens[[w]], function(z)paste0(x,":",y,";",z)))))))

      # get unique combinations per cell
      cells <- unique(ltemp[[cell]])
      cellcombos <- lapply(cells,function(x)
        unique(unlist(combos[ltemp[[cell]] == x])))
      #lcounts <- table(unlist(lapply(cellcombos,function(x)x)))
      lcounts <- table(unlist(cellcombos,function(x)x))
      max <- names(lcounts)[which.max(lcounts)]
      cvs <- unlist(lapply(combos,function(x)max %in% x))
      ltemp[cvs,][[subgroup]] <- lclone
      ltemp[cvs,]$vj_gene <- max
      
      # if a cell has the same combo for two rearrangements, only pick one
      # with the most ACTG characters
      rmseqs <- c()
      cell_counts <- table(ltemp[cvs,][[cell]])
      mcells <- names(cell_counts)[cell_counts > 1]
      for(cellname in mcells){
        ttemp <- dplyr::filter(ltemp,cvs & !!rlang::sym(cell) == cellname)
        ttemp$str_counts <-
          stringr::str_count(ttemp[[seq]],"[A|C|G|T]")
        # keep version with most non-N characters
        keepseq <- ttemp[[id]][which.max(ttemp$str_counts)]
        rmtemp <- ttemp[!ttemp[[id]] == keepseq,]
        rmseqs <- c(rmseqs,rmtemp[[id]])
      }
      include <- dplyr::filter(ltemp, cvs & !(!!rlang::sym(id) %in% rmseqs))
      leave <- dplyr::filter(ltemp,!cvs | (!!rlang::sym(id) %in% rmseqs))
      
      # find other cells still in ltemp and add as vj_alt_cell
      mcells <- unique(include[[cell]])
      for(cellname in mcells){
        if(cellname %in% leave[[cell]]){
          include[include[[cell]] == cellname,]$vj_alt_cell <-
            paste(paste0(leave[leave[[cell]] == cellname,][[v_call]],":",
                         leave[leave[[cell]] == cellname,][[j_call]]),
                  collapse=",")
        }
      }
      ld <- dplyr::bind_rows(ld,include)
      ltemp <- dplyr::filter(ltemp,!(!!rlang::sym(cell) %in% ltemp[cvs,][[!!cell]]))
      lclone <- lclone + 1
    }
    ld[[clone]] <- cloneid
    for(cellname in unique(hd_sc[[cell]])){
      if(cellname %in% ld[[cell]]){
        lclone <- ld[ld[[cell]] == cellname,][[subgroup]]
        hd_sc[hd_sc[[cell]] == cellname,][[subgroup]] <- lclone
        hd_sc[hd_sc[[cell]] == cellname,]$vj_gene <- ld[ld[[cell]] == cellname,]$vj_gene
        hd_sc[hd_sc[[cell]] == cellname,]$vj_alt_cell <-
          ld[ld[[cell]] == cellname,]$vj_alt_cell
      }
    }
    # now get the subgroup_id for the heavy chains lacking paired light chains
    comb <- dplyr::bind_rows(hd_sc,ld)
    comb[[subgroup]] <- as.integer(comb[[subgroup]])
    if(nrow(ld) != 0 & nrow(hd_bulk) != 0){
      for(sequence in 1:nrow(hd_bulk)){
        # CGJ 8/6/24 -- updated to do the padding on the temp sequences
        rating <- unlist(lapply(1:nrow(hd_sc), function(x){
          temp <- rbind(hd_sc[x,], hd_bulk[sequence,])
          if(nchar(temp[[seq]][1] != nchar(temp[[seq]][2]))){
            temp <- alakazam::padSeqEnds(temp[[seq]])
          }
          value <- alakazam::seqDist(temp[1], temp[2])
          return(value)
        }))
        rating <- as.numeric(rating)
        # row number of heavy chain only df with lowest seq dist
        proper_index <- which(rating == min(rating))
        if(length(proper_index) > 1){
          # find the subgroups that belong to the lowest seq dists
          subgroups <- hd_sc[[subgroup]][proper_index]
          if(length(unique(subgroups)) > 1){
            # if there is more than one subgroup find the subgroup sizes of the 
            # subgroups being considered
            subgroup_size <- data.frame(clone_subgroup = unique(subgroups))
            subgroup_size$sizes <- unlist(lapply(1:nrow(subgroup_size), function(x){
              nrow(hd_sc[hd_sc[[subgroup]] == subgroup_size$clone_subgroup[x],])
            }))
            # if there is one subgroup that is the largest use it
            if(length(which(subgroup_size$sizes == max(subgroup_size$sizes))) == 1){
              proper_index_value <- subgroup_size$clone_subgroup[
                which(subgroup_size$sizes == max(subgroup_size$sizes))]
            } else { 
              # if there are more than one subgroup with the same size use the lower number
              potential_subgroups <- subgroup_size$clone_subgroup[
                which(subgroup_size$sizes == max(subgroup_size$sizes))]
              proper_index_value <- min(potential_subgroups)
            }
          } else{
            proper_index_value <- hd_sc[[subgroup]][proper_index[1]]
          }
        } else{
          proper_index_value <- hd_sc[[subgroup]][proper_index]
        }
        hd_bulk[[subgroup]][sequence] <- proper_index_value
      }
    } 
    if(nrow(hd_bulk) != 0){
      comb <- dplyr::bind_rows(comb, hd_bulk)
    }
    comb$clone_subgroup_id <- paste0(comb[[clone]],"_",comb[[subgroup]])
    comb$vj_cell <- sapply(1:nrow(comb), function(x){
      if(!is.na(comb$vj_alt_cell[x])){
        paste(comb$vj_gene[x],comb$vj_alt_cell[x],sep=",")        
      }else{
        comb$vj_gene[x]
      }
    })
    
    size <- c()
    for(subgroups in sort(unique(comb[[subgroup]]))){
      size <- append(size, nrow(comb[comb[[subgroup]] == subgroups,]))
    }
    if(!all(diff(size) <= 0)){
      order_check <- data.frame(table(comb[[subgroup]]))
      colnames(order_check) <- c(subgroup, "size")
      order_check <- order_check[order(-order_check$size),]
      order_check$proper_subgroup <- 1:nrow(order_check)
      comb$new_subgroup <- NA
      for(i in unique(comb[[subgroup]])){
        comb$new_subgroup[comb[[subgroup]] == i] <- order_check$proper_subgroup[order_check[[subgroup]] == i]
      }
      comb <- comb[, -which(names(comb) == subgroup)]
      names(comb)[names(comb) == "new_subgroup"] <- subgroup
    }
    comb
  },mc.cores=nproc)
  paired <- dplyr::bind_rows(paired)
  # remove the junction length from the vj_gene
  paired$vj_gene <- gsub("\\;.", "", paired$vj_gene)
  return(paired)
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
    if(!inherits(clones$data[[1]], "airrClone")){
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
    x@data$sequence_id=gsub(":|;|,|=| ","_",x@data$sequence_id);
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
