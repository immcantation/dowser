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
#' @param    vcall        name of the column containing V-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    jcall        name of the column containing J-segment allele assignments. All 
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
#' @param     collapse    iollapse identical sequences?
#' @param     traits      column ids to keep distinct during sequence collapse 
#' @param     region      if HL, include light chain information if available.
#' @param     heavy       name of heavy chain locus (default = "IGH")
#' @param     cell        name of the column containing cell assignment information
#' @param     locus       name of the column containing locus information
#' @param     triple      pad sequences to length mutliple three?
#' @return   A \link{airrClone} object containing the modified clone.
#'
#' @details
#' The input data.frame (\code{data}) must columns for each of the required column name 
#' arguments: \code{id}, \code{seq}, \code{germ}, \code{vcall}, \code{jcall}, 
#' \code{junc_len}, and \code{clone}.  The default values are as follows:
#' \itemize{
#'   \item  \code{id       = "sequence_id"}:         unique sequence identifier.
#'   \item  \code{seq      = "sequence_alignment"}:  IMGT-gapped sample sequence.
#'   \item  \code{germ     = "germline_alignment"}:  IMGT-gapped germline sequence.
#'   \item  \code{vcall    = "v_call"}:              V segment allele call.
#'   \item  \code{jcall    = "j_call"}:              J segment allele call.
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
#' \code{germ}, \code{vcall}, \code{jcall}, \code{junc_len} and \code{clone} columns, 
#' respectively. For any given clone, each value in these columns should be identical.
#'  
#' @seealso  Executes in order \link{maskSeqGaps}, \link{maskSeqEnds}, 
#'           \link{padSeqEnds}, and \link{collapseDuplicates}. 
#'           Returns a \link{airrClone} object which serves as input to
#'           \link{buildPhylipLineage}.
#' 
#' @export
#requires one loci to be the "primary" which is present in all cells and 
#is assumed to descend from a single common ancestor via point mutations
#and allow for one other alternate loci, which is assumed to descend by
#point mutations from a common ancestor
makeAirrClone <- 
function(data, id="sequence_id", seq="sequence_alignment", 
    germ="germline_alignment", vcall="v_call", jcall="j_call",
    junc_len="junction_length", clone="clone_id", mask_char="N",
    max_mask=0, pad_end=FALSE, text_fields=NULL, num_fields=NULL, seq_fields=NULL,
    add_count=TRUE, verbose=FALSE, collapse=TRUE, region="H", heavy=NULL,
    cell="cell", locus="locus", traits=NULL, triple=FALSE, randomize=TRUE){

    # Check for valid fields
    check <- alakazam::checkColumns(data, unique(c(id, seq, germ, vcall, jcall, junc_len, clone, 
                                      text_fields, num_fields, seq_fields,traits)))
    if (check != TRUE) { stop(check) }

    if(region=="HL"){
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
        tmp_df <- data[, unique(c(id, seq, text_fields, num_fields, seq_fields, cell, locus, traits))]
        tmp_df[[seq]] <- alakazam::maskSeqGaps(tmp_df[[seq]], mask_char=mask_char, 
            outer_only=FALSE)
        hc <- dplyr::filter(tmp_df,!!rlang::sym(locus)==rlang::sym(heavy))
        alt <- dplyr::filter(tmp_df,!!rlang::sym(locus)!=rlang::sym(heavy))
        if(nrow(hc) == 0){
            stop(paste("clone",unique(dplyr::pull(data,clone)),
                "heavy chain locus not found in dataset!"))
        }
        if(nrow(alt) == 0){
            region <- "H"
        }else{
            if(length(unique(dplyr::pull(alt,!!locus))) > 1){
                stop(paste("clone",unique(dplyr::pull(data,clone)),
                    "currently only one alternate loci per clone supported"))
            }
        }
    }else{
        # Replace gaps with Ns and masked ragged ends
        tmp_df <- data[, unique(c(id, seq, text_fields, num_fields, seq_fields, traits))]
        tmp_df[[seq]] <- alakazam::maskSeqGaps(tmp_df[[seq]], mask_char=mask_char, 
            outer_only=FALSE)
    }

    if(region=="HL"){
        #print("Concatenating H and L")
        hc[[seq]] <- alakazam::maskSeqEnds(hc[[seq]], mask_char=mask_char, 
            max_mask=max_mask, trim=FALSE)
        alt[[seq]] <- alakazam::maskSeqEnds(alt[[seq]], mask_char=mask_char, 
            max_mask=max_mask, trim=FALSE)
        # Pad ends
        if(pad_end) {
            hc[[seq]] <- alakazam::padSeqEnds(hc[[seq]], pad_char=mask_char, triple=triple)
            alt[[seq]] <- alakazam::padSeqEnds(alt[[seq]], pad_char=mask_char, triple=triple)
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
        germline <- alakazam::maskSeqGaps(hcd[[germ]][1], mask_char=mask_char, outer_only=FALSE)
        lgermline <- alakazam::maskSeqGaps(altd[[germ]][1], mask_char=mask_char, outer_only=FALSE)
        if(pad_end){
             germline <- alakazam::padSeqEnds(germline, pad_char=mask_char, triple=triple)
            lgermline <- alakazam::padSeqEnds(lgermline, pad_char=mask_char, triple=triple)
        }
        hlgermline <- paste0(germline,lgermline)
        tmp_df <- hc
        loci <- unique(dplyr::pull(data,!!locus))
        tmp_df[[rlang::sym(locus)]] <- paste(loci,collapse=",")
        regions <- c(rep(unique(dplyr::pull(hc,!!locus)),times=hc_length),
                 rep(unique(dplyr::pull(alt,!!locus)),times=alt_length))
        numbers <- c(1:hc_length,1:alt_length)
        if(length(regions) != unique(nchar(tmp_df$hlsequence))){
            stop(paste("clone",unique(dplyr::pull(data,clone)),
                "regions vector not equal to total sequence length!"))
        }
        if(length(regions) != nchar(hlgermline)){
            stop(paste("clone",unique(dplyr::pull(data,clone)),
                "regions vector not equal to germline sequence length!"))
        }
        new_seq <- "hlsequence"
    }else{
        tmp_df[[seq]] <- alakazam::maskSeqEnds(tmp_df[[seq]], 
            mask_char=mask_char, max_mask=max_mask, trim=FALSE)
        if(pad_end){
            tmp_df[[seq]] <- alakazam::padSeqEnds(tmp_df[[seq]], pad_char=mask_char, triple=triple)
        }
        germline <- alakazam::maskSeqGaps(data[[germ]][1], 
            mask_char=mask_char, outer_only=FALSE)
        if(pad_end){
            germline <- alakazam::padSeqEnds(germline, pad_char=mask_char, triple=triple)
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
        regions <- rep(loci,times=nchar(germline))
        numbers <- 1:nchar(germline)
        lgermline <- ""
        hlgermline <- germline
        tmp_df$lsequence <- ""
        tmp_df$hlsequence <- tmp_df[[seq]]
        new_seq <- seq
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
        tmp_df = tmp_df[sample(1:nrow(tmp_df),replace=FALSE),]
    }
    
    # Define return object
    tmp_names <- names(tmp_df)
    if ("sequence" %in% tmp_names & seq != "sequence") {
        tmp_df <- tmp_df[, tmp_names != "sequence"]
        tmp_names <- names(tmp_df)
    }
    names(tmp_df)[tmp_names == seq] <- "sequence"
    names(tmp_df)[tmp_names == id] <- "sequence_id"
    
    if(region=="HL"){
        phylo_seq <- "hlsequence"
    }else if(region=="L"){
        phylo_seq <- "lsequence"
    }else{
        phylo_seq <- "sequence"
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
        v_gene=alakazam::getGene(data[[vcall]][1]), 
        j_gene=alakazam::getGene(data[[jcall]][1]), 
        junc_len=data[[junc_len]][1],
        locus=unique(loci),
        region=regions,
        numbers=numbers,
        phylo_seq=phylo_seq)
    
    outclone
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
#' @param    data         data.frame containing the AIRR or Change-O data for a clone. See Details
#'                        for the list of required columns and their default values.
#' @param    id           name of the column containing sequence identifiers.
#' @param    seq          name of the column containing observed DNA sequences. All 
#'                        sequences in this column must be multiple aligned.
#' @param    subclones    split or lump subclones? See \code{getSubclones}.
#' @param    minseq       minimum numbner of sequences per clone
#' @param    triple       pad sequences to length multiple of three
#' @param    majoronly    only return largest subclone and sequences without light chains
#' @param    germ         name of the column containing germline DNA sequences. All entries 
#'                        in this column should be identical for any given clone, and they
#'                        must be multiple aligned with the data in the \code{seq} column.
#' @param    vcall        name of the column containing V-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    jcall        name of the column containing J-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    junc_len     name of the column containing the length of the junction as a 
#'                        numeric value. All entries in this column should be identical 
#'                        for any given clone.
#' @param    clone        name of the column containing the identifier for the clone. All 
#'                        entries in this column should be identical.
#' @param    subclone     name of the column containing the identifier for the subclone.
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
#' @param     collapse    iollapse identical sequences?
#' @param     traits      column ids to keep distinct during sequence collapse 
#' @param     region      if HL, include light chain information if available.
#' @param     heavy       name of heavy chain locus (default = "IGH")
#' @param     cell        name of the column containing cell assignment information
#' @param     locus       name of the column containing locus information
#' @param     nproc		  Number of cores to parallelize formating over.                        
#'
#' @return   A list of \link{airrClone} objects containing modified clones.
#'
#' @details
#' This function is largely a wrapper for alakazam::makeAirrClone.
#' The input data.frame (\code{data}) must columns for each of the required column name 
#' arguments: \code{id}, \code{seq}, \code{germ}, \code{vcall}, \code{jcall}, 
#' \code{junc_len}, and \code{clone}.  The default values are as follows:
#' \itemize{
#'   \item  \code{id       = "sequence_id"}:           unique sequence identifier.
#'   \item  \code{seq      = "sequence_alignment"}:         IMGT-gapped sample sequence.
#'   \item  \code{germ     = "germline_alignment_d_mask"}:  IMGT-gapped germline sequence.
#'   \item  \code{vcall    = "v_call"}:                V-segment allele call.
#'   \item  \code{jcall    = "j_call"}:                J-segment allele call.
#'   \item  \code{junc_len = "junction_length"}:       junction sequence length.
#'   \item  \code{clone    = "clone_id"}:                 clone identifier.
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
#' \code{germ}, \code{vcall}, \code{jcall}, \code{junc_len} and \code{clone} columns, 
#' respectively. For any given clone, each value in these columns should be identical.
#'  
#' @seealso  Executes in order \code{alakazam::maskSeqGaps}, \code{alakazam::maskSeqEnds}, 
#'           \code{alakazam::padSeqEnds}, \code{alakazam::collapseDuplicates},
#' 			 and \code{processClones}. Returns a list of \link{airrClone} objects 
#' 			which serve as input to \link{getTrees} and \link{bootstrapTrees}.
#' 
#' @examples
#' \dontrun{
#' data(ExampleDb)
#' clones <- formatClones(ExampleDb,trait="sample_id")
#' }
#' @export
formatClones <- function(data, id="sequence_id", seq="sequence_alignment", 
                germ="germline_alignment_d_mask", vcall="v_call", jcall="j_call",
                junc_len="junction_length", clone="clone_id", subclone="subclone_id",
                mask_char="N", max_mask=0, pad_end=TRUE, text_fields=NULL, num_fields=NULL, 
                seq_fields=NULL, add_count=TRUE, verbose=FALSE, nproc=1, collapse=TRUE,
                region="H", heavy=NULL, cell="cell_id", locus="locus", traits=NULL,
                minseq=2, triple=TRUE, subclones="lump", majoronly=FALSE,randomize=TRUE) {

	if(majoronly){
		if(!subclone %in% names(data)){
			stop("Need subclone designation if majoronly=TRUE")
		}
		data <- filter(data, !!rlang::sym(subclone) <= 1)
	}

	if(region == "H"){ #if region is heavy and, discard all non-IGH sequences
		if(!is.null(heavy)){
			if(locus %in% names(data)){
				data <- filter(data, !!rlang::sym(locus) == rlang::sym(heavy))
			}
		}
	}
	if(region == "HL"){
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
	if(subclones == "split"){
		if(!subclone %in% names(data)){
			stop("Need subclone designation for heavy+light chain clones")
		}
		if(sum(data[[subclone]] == 0) > 0){
			warning("Assigning subclone 0 (missing light chain) to subclone 1")
			data[data[[subclone]] == 0,][[subclone]] <- 1
		}
		data[[clone]] <- paste0(data[[clone]],"_",data[[subclone]])
	}else if(subclones == "lump" && region=="HL"){
		data <- filter(data, !(!!rlang::sym(locus) != rlang::sym(heavy) &
			!!rlang::sym(subclone) > 1))
	}else if(region == "HL"){
		stop("subclones designation must be either lump or split")
	}

    clones <- data %>%
        dplyr::group_by(!!rlang::sym(clone)) %>%
        dplyr::do(data=makeAirrClone(.,
            id=id, seq=seq, germ=germ, vcall=vcall, jcall=jcall, junc_len=junc_len,
            clone=clone, mask_char=mask_char, max_mask=max_mask, pad_end=pad_end,
            text_fields=text_fields, num_fields=num_fields, seq_fields=seq_fields,
            add_count=add_count, verbose=verbose,collapse=collapse,
            region=region, heavy=heavy, cell=cell, locus=locus, traits=traits,triple=triple,
            randomize=randomize))

    if(region == "HL"){
    	seq_name <- "hlsequence"
    }else{
    	seq_name <- "sequence"
    }
    
    fclones <- processClones(clones, nproc=nproc, seq=seq_name, minseq=minseq)
    fclones
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
#' @param    cellid       name of the column containing identifier for cells.
#' @param    vcall        name of the column containing V-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    jcall        name of the column containing J-segment allele assignments. All 
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
#' TODO: Junction length?
#' TODO: Option to just store all VJ pairs for a cell in the heavy, remove light seqs
#' TODO: Option to split VJ parititions into separate clones
#' TODO: Make v_alt_cell not be NA by default
#' TODO: Unit testing
#' @export
getSubclones <- function(heavy, light, nproc=1, minseq=2,
	id="sequence_id", seq="sequence_alignment", 
	clone="clone_id", cellid="cell_id", vcall="v_call", jcall="j_call",
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
		#print(cloneid)
		hd <- filter(heavy,!!rlang::sym(clone) == cloneid)
		ld <- filter(light,!!rlang::sym(cellid) %in% hd[[cellid]])
		hd <- filter(hd,(!!rlang::sym(cellid) %in% ld[[cellid]]))
		hr <- filter(hd,!(!!rlang::sym(cellid) %in% ld[[cellid]]))
		if(nrow(ld) == 0){
			return(hd)
		}
		ltemp <- ld
		ltemp$clone_id <- -1
		ld <- dplyr::tibble()
		lclone <- 1
		while(nrow(ltemp) > 0){
			lvs <- strsplit(ltemp[[vcall]],split=",")
			ljs <- strsplit(ltemp[[jcall]],split=",")
			combos <- 
				lapply(1:length(lvs),function(w)
				unlist(lapply(lvs[[w]],function(x)
				lapply(ljs[[w]],function(y)paste(x,y,sep=":")))))
			cells <- unique(ltemp[[cellid]])
			cellcombos <- lapply(cells,function(x)
				unique(unlist(combos[ltemp[[cellid]] == x])))
			lcounts <- table(unlist(lapply(cellcombos,function(x)x)))
			max <- names(lcounts)[which.max(lcounts)]
			cvs <- unlist(lapply(combos,function(x)max %in% x))
			ltemp[cvs,][[subclone]] <- lclone
			ltemp[cvs,]$vj_gene <- max

			# if a cell has the same combo for two rearrangements, only pick one
			rmseqs <- c()
			cell_counts <- table(ltemp[cvs,][[cellid]])
			mcells <- names(cell_counts)[cell_counts > 1]
			for(cell in mcells){
				ttemp <- filter(ltemp,cvs & !!rlang::sym(cellid) == cell)
				ttemp$str_counts <- 
				    stringr::str_count(ttemp[[seq]],"[A|C|G|T]")
				rmtemp <- ttemp[-which.max(ttemp$str_counts),]
				rmseqs <- c(rmseqs,rmtemp[[id]])
			}
			include <- filter(ltemp,cvs & !(!!rlang::sym(id) %in% rmseqs))
			leave <- filter(ltemp,!cvs | (!!rlang::sym(id) %in% rmseqs))

			# find other cells still in ltemp and add as vj_alt_cell
			mcells <- unique(include[[cellid]])
			for(cell in mcells){
				if(cell %in% leave[[cellid]]){
					include[include[[cellid]] == cell,]$vj_alt_cell <- 
						paste(paste0(leave[leave[[cellid]] == cell,][[vcall]],":",
							leave[leave[[cellid]] == cell,][[jcall]]),
							collapse=",")
				}
			}
			ld <- bind_rows(ld,include)
			ltemp <- filter(ltemp,!(!!rlang::sym(cellid) %in% ltemp[cvs,][[cellid]]))
			lclone <- lclone + 1
		}
		ld[[clone]] <- cloneid
		for(cell in unique(hd[[cellid]])){
			hclone <- hd[hd[[cellid]] == cell,][[clone]]
			if(cell %in% ld[[cellid]]){
				lclone <- ld[ld[[cellid]] == cell,][[subclone]]
				ld[ld[[cellid]] == cell,][[subclone]] <- lclone
				hd[hd[[cellid]] == cell,][[subclone]] <- lclone
				hd[hd[[cellid]] == cell,]$vj_gene <- ld[ld[[cellid]] == cell,]$vj_gene
				hd[hd[[cellid]] == cell,]$vj_alt_cell <- ld[ld[[cellid]] == cell,]$vj_alt_cell
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

#' Wrapper for CreateGermlines.py
#' 
#' \code{createGermlines} reconstructs clonal and subclonal germlines
#' 
#' @param    data         tibble containing sequence information
#' @param    exec         location of CreateGermline.py
#' @param    refs         vector of reference allele locations
#' @param    file         temporary file name to write
#' @param    cf           column name for clone or subclone id
#' @param    format       airr or changeo format
#' @param    g            full or dmask germline option
#' @param    germ         name of the column containg germline DNA sequences in
#'                        output file. 
#' @param    rm_file      remove temporary file?
#'
#' @return   a tibble containing reconstructed germlines
#'  
#' @export
createGermlines <- function(data, exec, refs, file, cf="vj_clone",
		format="airr",g="dmask",germ="germline_alignment_d_mask",rm_file=TRUE){
	alakazam::writeChangeoDb(data,file=paste0(file,".tsv"))
	exec <- path.expand(exec)
	r <- paste(path.expand(refs),collapse=" ")
	command <- paste("-d",paste0(file,".tsv"),"-r",r,"--format",format,"--cloned --cf",cf,
		"--outname",file,"-g",g)
	params <- list(exec,command,stdout=TRUE,stderr=TRUE)
	status <- tryCatch(do.call(base::system2, params), error=function(e){
		return(e)
		}, warning=function(w){
		return(w)
		})
	if(class(status) != "character"){
		print(paste(exec,command))
		stop(status)
	}
	gl <- alakazam::readChangeoDb(paste0(file,"_germ-pass.tsv"))
	if(rm_file){
		unlink(file)
		unlink(paste0(file,"_germ-pass.tsv"))
	}
	if(g != "dmask"){
		if("germline_alignment" %in% names(gl)){
			gl[[germ]] <- gl$germline_alignment
		}else if("GERMLINE_IMGT" %in% names(gl)){
			gl[[germ]] <- gl$GERMLINE_IMGT
		}else{
			stop("New germline column not found in input file")
		}
	}
	return(gl)
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
		x@data$sequence_id=	gsub(":","_",x@data$sequence_id);
		x })
	clones$data <- lapply(clones$data,function(x){
		x@data$sequence_id=gsub(";","_",x@data$sequence_id);
		x })
	
	clones$data <- lapply(clones$data,function(x){
		x@data$sequence_id=gsub(",","_",x@data$sequence_id);
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
		clones$locus <- unlist(lapply(clones$data,function(x)paste(x@locus,collapse=",")))
	}
	clones$seqs <- unlist(lapply(clones$data,function(x)nrow(x@data)))
	clones <- dplyr::rowwise(clones)
	clones <- dplyr::ungroup(clones)
	clones
}