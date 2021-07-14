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
#' @param    chain        if HL, include light chain information if available.
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
#' \code{junc_len}, and \code{clone}.  The default values are as follows:
#' \itemize{
#'   \item  \code{id       = "sequence_id"}:         unique sequence identifier.
#'   \item  \code{seq      = "sequence_alignment"}:  IMGT-gapped sample sequence.
#'   \item  \code{germ     = "germline_alignment"}:  IMGT-gapped germline sequence.
#'   \item  \code{v_call    = "v_call"}:              V segment allele call.
#'   \item  \code{j_call    = "j_call"}:              J segment allele call.
#'   \item  \code{junc_len = "junction_length"}:     junction sequence length.
#'   \item  \code{clone    = "clone_id"}:            clone identifier.
#' }
#' Additional annotation columns specified in the \code{text_fields}, \code{num_fields} 
#' or \code{seq_fields} arguments will be retained in the \code{data} slot of the return 
#' object, but are not required. If the input data.frame \code{data} already contains a 
#' column named \code{sequence}, which is not used as the \code{seq} argument, then that 
#' column will not be retained.
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
#' @seealso  Returns an \link{airrClone}. See \link{formatClones} to enerate an 
#' ordered list of airrClone objects.
#' @example
#' airr_clone <- makeAirrClone(ExampleDb[ExampleDb$clone_id=="3184",])
#' @export
#requires one loci to be the "primary" which is present in all cells and 
#is assumed to descend from a single common ancestor via point mutations
#and allow for one other alternate loci, which is assumed to descend by
#point mutations from a common ancestor
makeAirrClone <- 
function(data, id="sequence_id", seq="sequence_alignment", 
    germ="germline_alignment_d_mask", v_call="v_call", j_call="j_call",
    junc_len="junction_length", clone="clone_id", mask_char="N",
    max_mask=0, pad_end=TRUE, text_fields=NULL, num_fields=NULL, seq_fields=NULL,
    add_count=TRUE, verbose=FALSE, collapse=TRUE, chain="H", heavy=NULL,
    cell="cell_id", locus="locus", traits=NULL, mod3=TRUE, randomize=TRUE,
    use_regions=TRUE, dup_singles=FALSE){

    # Check for valid fields
    check <- alakazam::checkColumns(data, 
        unique(c(id, seq, germ, v_call, j_call, junc_len, clone, 
        text_fields, num_fields, seq_fields,traits)))
    if (check != TRUE) { stop(check) }

    if(chain=="HL"){
        check <- alakazam::checkColumns(data, c(cell,locus))
        if (check != TRUE) { stop(check) }
        if(is.null(heavy)){
            stop(paste("clone",unique(dplyr::pull(data,clone)),
                "heavy chain loci ID must be specified if combining loci!"))
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
        if(length(hc_length) > 1){
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
        }else{
            regions <- rep("N", times=nchar(hlgermline))
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
        numbers <- 1:nchar(germline)
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
        if (!pad_end) {
            len_message <- paste(len_message, 
                "Consider specifying pad_end=TRUE and verify the multiple alignment.")
        } else {
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
        clone=as.character(data[[clone]][1]),
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
        #chain=chains,
        region=regions,
        numbers=numbers,
        phylo_seq=phylo_seq)
    
    outclone
}


# Remove uniformative columns from data and germline
# 
# \code{cleanAlignment} clean multiple sequence alignments
# @param    clone   \code{airrClone} object
# @param    seq     column in \code{clone} object
#
# @return   \code{airrClone} object with cleaned alignment
#
cleanAlignment <- function(clone, seq="sequence"){
    if(seq=="hlsequence"){
        g <- strsplit(clone@hlgermline[1],split="")[[1]]
    }else{
        g <- strsplit(clone@germline[1],split="")[[1]]
    }
    sk <- strsplit(clone@data[[seq]],split="")
    sites=seq(1,length(g)-3,by=3)
    ns <- c()
    for(i in sites){
        l=lapply(sk,function(x) paste(x[i:(i+2)],collapse="")=="NNN")
        ns <- c(ns,sum(unlist(l)),sum(unlist(l)),sum(unlist(l)))
    }
    informative <- ns != length(sk)
    l=lapply(sk,function(x) x=paste(x[informative],collapse=""))
    gm=paste(g[informative],collapse="")
    if(.hasSlot(clone,"locus")){
        clone@locus <- clone@locus[informative] 
    }
    if(.hasSlot(clone,"region")){
        clone@region <- clone@region[informative]   
    }
    if(.hasSlot(clone,"numbers")){
        clone@numbers <- clone@numbers[informative]
    }
    if(seq=="hlsequence"){
        clone@hlgermline=gm
    }else{
        clone@germline=gm
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
#' @param    chain        if HL, include light chain information if available.
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
#'      which serve as input to \link{getTrees} and \link{bootstrapTrees}.
#' 
#' @examples
#' data(ExampleDb)
#' # Select two clones, for demonstration purpose
#' sel <- c("3170", "3184")
#' clones <- formatClones(ExampleDb[ExampleDb$clone_id %in% sel,],trait="sample_id")
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
    }

    counts <- table(data[[clone]])
    rmclones <- names(counts[counts < minseq])
    data <- data[!data[[clone]] %in% rmclones,]

    data <- data[!is.na(data[[seq]]),]
    clones <- data %>%
        dplyr::group_by(!!rlang::sym(clone)) %>%
        dplyr::do(data=makeAirrClone(.data, seq=seq,
            clone=clone, chain=chain, heavy=heavy, cell=cell,...))

    if(chain == "HL"){
        seq_name <- "hlsequence"
    }else{
        seq_name <- "sequence"
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


#' \code{maskCodons} Masks codons split by insertions
#' @param    id              sequence id
#' @param    q               (query) un-aligned input sequence (sequence)
#' @param    s               (subject) aligned input sequence (sequence_alignment)
#' @param    keep_alignment  store q and s alignments
#' @param    keep_insertions return removed insertion sequences?
#' @param    gap_opening      gap opening penalty (Biostrings::pairwiseALignment)
#' @param    gap_extension    gap extension penalty (Biostrings::pairwiseALignment)
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


#' Define subclones based on light chain rearrangements
#' 
#' \code{getSubclones} plots a tree or group of trees
#' @param    heavy        A tibble containing heavy chain sequences with clone_id
#' @param    light        A tibble containing light chain sequences
#' @param    nproc        number of cores for parallelization
#' @param    minseq       minimum number of sequneces per clone
#' @param    id           name of the column containing sequence identifiers.
#' @param    seq          name of the column containing observed DNA sequences. All 
#'                        sequences in this column must be multiple aligned.
#' @param    clone        name of the column containing the identifier for the clone. All 
#'                        entries in this column should be identical.
#' @param    cell_id       name of the column containing identifier for cells.
#' @param    v_call        name of the column containing V-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    j_call        name of the column containing J-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    junc_len     name of the column containing the length of the junction as a 
#'                        numeric value. All entries in this column should be identical 
#'                        for any given clone.
#' @param    nolight      string to use to indicate a missing light chain
#'
#' @return   a tibble containing 
#'
#' @details
#' 1. Make temporary array containing light chain clones
#' 2. Enumerate all possible V and J combinations
#' 3. Determine which combination is the most frequent
#' 4. Assign sequences with that combination to clone t
#' 5. Copy those sequences to return array
#' 6. Remove all cells with that combination from temp array
#' 7. Repeat 1-5 until temporary array zero.
#' If there is more than rearrangement with the same V/J
#' in the same cell, pick the one with the highest non-ambiguous
#' characters. 
# TODO: Junction length?
# TODO: Option to just store all VJ pairs for a cell in the heavy, remove light seqs
# TODO: Option to split VJ parititions into separate clones
# TODO: Make v_alt_cell not be NA by default
#' @export
getSubclones <- function(heavy, light, nproc=1, minseq=1,
    id="sequence_id", seq="sequence_alignment", 
    clone="clone_id", cell_id="cell_id", v_call="v_call", j_call="j_call",
    junc_len="junction_length", nolight="missing"){
    
    subclone <- "subclone_id"
    scount <- table(heavy[[clone]])
    big <- names(scount)[scount >= minseq]
    heavy <- filter(heavy,(!!rlang::sym(clone) %in% big))

    heavy$vj_gene <- nolight
    heavy$vj_alt_cell <- nolight
    heavy$subclone_id <- 0
    light$vj_gene <- nolight
    light$vj_alt_cell <- nolight
    light$subclone_id <- 0
    light[[clone]] <- -1
    paired <- parallel::mclapply(unique(heavy[[clone]]),function(cloneid){
        hd <- filter(heavy,!!rlang::sym(clone) == cloneid)
        ld <- filter(light,!!rlang::sym(cell_id) %in% hd[[!!cell_id]])
        hd <- filter(hd,(!!rlang::sym(cell_id) %in% ld[[!!cell_id]]))
        # hr <- filter(hd,!(!!rlang::sym(cell_id) %in% ld[[!!cell_id]]))
        if(nrow(ld) == 0){
            return(hd)
        }
        ltemp <- ld
        ltemp$clone_id <- -1
        ld <- dplyr::tibble()
        lclone <- 1
        while(nrow(ltemp) > 0){
            lvs <- strsplit(ltemp[[v_call]],split=",")
            ljs <- strsplit(ltemp[[j_call]],split=",")
            combos <- 
                lapply(1:length(lvs),function(w)
                unlist(lapply(lvs[[w]],function(x)
                lapply(ljs[[w]],function(y)paste(x,y,sep=":")))))
            cells <- unique(ltemp[[cell_id]])
            cellcombos <- lapply(cells,function(x)
                unique(unlist(combos[ltemp[[cell_id]] == x])))
            lcounts <- table(unlist(lapply(cellcombos,function(x)x)))
            max <- names(lcounts)[which.max(lcounts)]
            cvs <- unlist(lapply(combos,function(x)max %in% x))
            ltemp[cvs,][[subclone]] <- lclone
            ltemp[cvs,]$vj_gene <- max

            # if a cell has the same combo for two rearrangements, only pick one
            rmseqs <- c()
            cell_counts <- table(ltemp[cvs,][[cell_id]])
            mcells <- names(cell_counts)[cell_counts > 1]
            for(cell in mcells){
                ttemp <- filter(ltemp,cvs & !!rlang::sym(cell_id) == cell)
                ttemp$str_counts <- 
                    stringr::str_count(ttemp[[seq]],"[A|C|G|T]")
                rmtemp <- ttemp[-which.max(ttemp$str_counts),]
                rmseqs <- c(rmseqs,rmtemp[[id]])
            }
            include <- filter(ltemp,cvs & !(!!rlang::sym(id) %in% rmseqs))
            leave <- filter(ltemp,!cvs | (!!rlang::sym(id) %in% rmseqs))

            # find other cells still in ltemp and add as vj_alt_cell
            mcells <- unique(include[[cell_id]])
            for(cell in mcells){
                if(cell %in% leave[[cell_id]]){
                    include[include[[cell_id]] == cell,]$vj_alt_cell <- 
                        paste(paste0(leave[leave[[cell_id]] == cell,][[v_call]],":",
                            leave[leave[[cell_id]] == cell,][[j_call]]),
                            collapse=",")
                }
            }
            ld <- bind_rows(ld,include)
            ltemp <- filter(ltemp,!(!!rlang::sym(cell_id) %in% ltemp[cvs,][[!!cell_id]]))
            lclone <- lclone + 1
        }
        ld[[clone]] <- cloneid
        for(cell in unique(hd[[cell_id]])){
            #hclone <- hd[hd[[cell_id]] == cell,][[clone]]
            if(cell %in% ld[[cell_id]]){
                lclone <- ld[ld[[cell_id]] == cell,][[subclone]]
                ld[ld[[cell_id]] == cell,][[subclone]] <- lclone
                hd[hd[[cell_id]] == cell,][[subclone]] <- lclone
                hd[hd[[cell_id]] == cell,]$vj_gene <- ld[ld[[cell_id]] == cell,]$vj_gene
                hd[hd[[cell_id]] == cell,]$vj_alt_cell <- ld[ld[[cell_id]] == cell,]$vj_alt_cell
            }
        }
        comb <- bind_rows(hd,ld)
        comb$vj_clone <- paste0(comb[[clone]],"_",comb[[subclone]])
        comb$vj_cell <- paste(comb$vj_gene,comb$vj_alt_cell,sep=",")
        comb
    },mc.cores=nproc)
    paired <- bind_rows(paired)
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

    or <- order(unlist(lapply(clones$data,function(x) nrow(x@data))),
        decreasing=TRUE)
    clones <- clones[or,]

    clones$data <- parallel::mclapply(clones$data,
        function(x)cleanAlignment(x,seq),mc.cores=nproc)

    if(.hasSlot(clones$data[[1]],"locus")){
        clones$locus <- unlist(lapply(clones$data,function(x)
            paste(sort(unique(x@locus)),collapse=",")))
    }
    clones$seqs <- unlist(lapply(clones$data,function(x)nrow(x@data)))
    clones <- dplyr::rowwise(clones)
    clones <- dplyr::ungroup(clones)
    clones
}