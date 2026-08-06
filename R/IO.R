# Input/output related functions

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


# Write a clone's sequence alignment to a fasta file
# 
# \code{cloneToFasta} write clone sequences as a fasta file
# @param    c            airrClone object
# @param    fastafile    file to be exported
# @param    germid       sequence id of germline
# @param    trait        trait to include in sequence ids
# @param    empty        don't include real sequence information
#
# @return   Name of exported fasta file.
cloneToFasta <- function(c, fastafile, germid, trait=NULL, empty=FALSE){
  text <- ""
  if(!is.null(trait)){
    c@data$sequence_id <- paste(c@data$sequence_id,c@data[,trait],sep="_")
  }
  for(i in 1:nrow(c@data)){
    text <- paste0(text,">",c@data[i,]$sequence_id,"\n")
    if(!empty){
      if(c@phylo_seq == "sequence"){
        text <- paste0(text,c@data[i,]$sequence,"\n")
      }else if(c@phylo_seq == "lsequence"){
        text <- paste0(text,c@data[i,]$lsequence,"\n")
      }else if(c@phylo_seq == "hlsequence"){
        text <- paste0(text,c@data[i,]$hlsequence,"\n")
      }else{
        stop(paste("phylo_seq not recognized",c@clone))
      }
    }else{
      text <- paste0(text,"ATG\n")
    }
  }
  text <- paste0(text,">",germid,"\n")
  if(!empty){
    if(c@phylo_seq == "sequence"){
      text <- paste0(text,c@germline,"\n")
    }else if(c@phylo_seq == "lsequence"){
      text <- paste0(text,c@lgermline,"\n")
    }else if(c@phylo_seq == "hlsequence"){
      text <- paste0(text,c@hlgermline,"\n")
    }else{
      stop(paste("phylo_seq not recognized",c@clone))
    }
  }else{
    text <- paste0(text,"ATG\n")
  }
  write(text,file=fastafile,append=FALSE)
  return(fastafile)
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


#' Read a fasta file into a list of sequences
#' \code{readFasta} reads a fasta file
#' @param    file      FASTA file
#'
#' @return   List of sequences
#' @export
readFasta <- function(file){
  f <- readLines(file)
  if(length(f) == 1){
    return(NULL)
  }
  seqs <- list()
  id <- NA
  for(line in f){
    if(grepl("^>",line)){
      id <- gsub(">","",line)
      seqs[[id]] <- ""
    }else{
      if(is.na(id)){
        stop(paste("Error reading",file))
      }
      seqs[[id]] <- paste0(seqs[[id]],line)
    }
  }
  seqs
}

#' Exports the phylogenetic trees from the airrClone object
#' 
#' \code{exportTrees}   Exports phylogenetic trees
#' @param clones         tibble \code{airrClone} objects, the output of 
#'                      \link{formatClones}
#' @param filepath      The file path for where the trees will be saved
#' @param tree_column   The name of the column that contains the trees
#' @param ...           additional arguments to be passed
#'  
#' @export
exportTrees <- function(clones, filepath, tree_column = "trees", ...){
  # check to see if the trees column is there 
  if(alakazam::checkColumns(clones, tree_column)){
    ape::write.tree(phy = clones$trees, file = filepath, ...)
  } else{
    stop(paste(tree_column, "not found in the input airrClone object. Please",
               "specify what column contains the phylogenetic trees."))
  }
  
}

#' Write the sequences used in tree building to a fasta format. If there are more 
#' than one tree in airrClone output the sequence id will be followed by "|clone_id".
#' 
#' \code{writeCloneSequences}   Exports the sequences used in tree building. 
#' @param clones         tibble \code{airrClone} objects, the output of 
#'                      \link{formatClones}
#' @param file       The file path and name of where the sequences will be saved
#'  
#' @export
writeCloneSequences <- function(clones, file){
  for(i in 1:nrow(clones)){
    clone_id <- clones$clone_id[i]
    # grab the germline 
    if(clones$data[[i]]@phylo_seq == "sequence"){
      germline <- clones$data[[i]]@germline
    } else if(clones$data[[i]]@phylo_seq == "hlsequence"){
      germline <- clones$data[[i]]@hlgermline
    } else if(clones$data[[i]]@phylo_seq == "lsequence"){
      germline <- clones$data[[i]]@lgermline
    }
    write(paste0(">Germline|", clone_id), file, append = TRUE)
    write(germline, file, append = TRUE)
    for(j in 1:nrow(clones$data[[i]]@data)){
      seq_id <- clones$data[[i]]@data$sequence_id[j]
      if(clones$data[[i]]@phylo_seq == "sequence"){
        sequence <- clones$data[[i]]@data$sequence[j]
      } else if(clones$data[[i]]@phylo_seq == "hlsequence"){
        sequence <- clones$data[[i]]@data$hlsequence[j]
      } else if(clones$data[[i]]@phylo_seq == "lsequence"){
        sequence <- clones$data[[i]]@data$lsequence[j]
      }
      write(paste0(">", seq_id, "|", clone_id), file, append = TRUE)
      write(sequence, file, append = TRUE)
    }
  }
}


#'\code{pmlParamsToList}
#' Convert a \code{phangorn::pml} model fit into a flat, JSON-serializable
#' list of fitted model parameters.
#' @param    fit  object of class \code{pml}, as stored in \code{tree$parameters}
#'                after \code{getTrees(..., build="pml")}
#' @noRd
pmlParamsToList <- function(fit){
  list(
    model     = fit$model,
    logLik    = fit$logLik,
    df        = fit$df,
    k         = fit$k,
    shape     = fit$shape,
    inv       = fit$inv,
    rate      = fit$rate,
    site.rate = fit$site.rate,
    bf        = stats::setNames(as.numeric(fit$bf), c("a","c","g","t")),
    Q         = stats::setNames(as.numeric(fit$Q), c("a_c","a_g","a_t","c_g","c_t","g_t"))
  )
}

#'\code{listToPmlParams}
#' Inverse of \code{pmlParamsToList}: rebuild a full, working
#' \code{phangorn::pml} object from a flat parameter list (as reconstructed
#' from a JSON file written by \code{writeTreesJSON}), the clone's real tree,
#' and its real tip sequences.
#' @param    lst   flat parameter list, as produced from \code{pmlParamsToList}
#' @param    tree  the clone's reconstructed tree with tip labels matching
#'                 \code{names(seqs)}
#' @param    seqs  named character vector of each tip's ungapped
#'                 sequence (including the germline tip), names matching
#'                 \code{tree$tip.label}
#' @noRd
listToPmlParams <- function(lst, tree, seqs){
  split_seqs <- lapply(seqs, function(s) strsplit(s, "")[[1]])
  if(length(lst$bf) == 20){
    real_data <- phangorn::phyDat(
      ape::as.AAbin(t(as.matrix(dplyr::bind_rows(split_seqs)))),
      type="AA")
  }else{
    real_data <- phangorn::phyDat(
      ape::as.DNAbin(t(as.matrix(dplyr::bind_rows(split_seqs)))))
  }

  fit <- phangorn::pml(
    tree      = ape::unroot(tree),
    data      = real_data,
    bf        = lst$bf,
    Q         = lst$Q,
    inv       = lst$inv,
    k         = lst$k,
    shape     = lst$shape,
    rate      = lst$rate,
    model     = lst$model,
    site.rate = lst$site.rate
  )
  # pml() recomputes logLik from bf/Q, which might be close but not identical
  # to the value stored at write time.
  fit$logLik <- lst$logLik
  fit
}

#'\code{dataFrameToRecords}
#' Convert a data.frame to a list of row-records (one named list per row)
#' for JSON serialization, converting any non-finite numeric value
#' (\code{NaN}, \code{Inf}, \code{-Inf}) to its string form first.
#'
#' \code{jsonlite::write_json(..., auto_unbox=TRUE)} silently drops the key
#' entirely for a data.frame cell holding \code{NaN} when the data.frame is
#' passed to it directly (not rendered as \code{null} -- the field is simply
#' absent from that row's JSON object). BEAST parameter logs routinely
#' contain genuine \code{NaN} values (e.g. \code{geometric.mean} for
#' non-log-scale parameters), so this must be handled explicitly rather than
#' relying on \code{write_json}'s default data.frame conversion.
#' Pre-converting non-finite values to their string form (each field stays a
#' length-1 scalar, so \code{auto_unbox} still unboxes it normally) avoids
#' this; \code{as.numeric()} on read correctly restores
#' \code{NaN}/\code{Inf}/\code{-Inf} from their string forms.
#' @param    df  a data.frame
#' @noRd
dataFrameToRecords <- function(df){
  lapply(seq_len(nrow(df)), function(i){
    row <- as.list(df[i, , drop=FALSE])
    lapply(row, function(v) if(is.numeric(v) && !is.finite(v)) as.character(v) else v)
  })
}


#'\code{writePosteriorTreeList}
#' Store one posterior tree list (\code{tr@info$trees_posterior} or
#' \code{tr@info$trees_with_traits_posterior}) into \code{clone$info}, under
#' field names built from \code{prefix} -- the same raw-treetext + tip-label
#' mechanism \code{\link{writeTimeTreesJSON}} uses for the single summary
#' tree, applied once per posterior sample. Used identically for both
#' posterior tree fields, since they're structurally the same kind of object
#' (a list of \code{treedata} samples) and may be present independently of
#' each other.
#' @param    clone_info  the in-progress \code{clone$info} list to add to
#' @param    posterior_trees  a list of \code{treedata} objects (one
#'                            \code{tr@info$trees_posterior} or
#'                            \code{tr@info$trees_with_traits_posterior})
#' @param    prefix      field name prefix, e.g. \code{"trees_posterior"} --
#'                       produces \code{clone_info[[prefix]]} (treetext list)
#'                       plus either \code{clone_info[[paste0(prefix,
#'                       "_tip_labels_shared")]]} or \code{[[paste0(prefix,
#'                       "_tip_labels")]]}
#' @param    ancestor_id           \code{"Germline"}
#' @param    ancestor_id_qualified the clone-qualified germline tip name
#' @param    has_germline          does the tree have a germline tip?
#' @noRd
writePosteriorTreeList <- function(clone_info, posterior_trees, prefix,
  ancestor_id, ancestor_id_qualified, has_germline){
  clone_info[[prefix]] <- lapply(posterior_trees, function(t) t@treetext)
  tip_label_list <- lapply(posterior_trees, function(t){
    labs <- t@phylo$tip.label
    if(has_germline){
      labs[labs == ancestor_id] <- ancestor_id_qualified
    }
    labs
  })
  # One shared tip-label array for all samples when they share one order
  # (the norm -- confirmed on real objects), else one per sample.
  shared <- length(tip_label_list) <= 1 ||
    all(vapply(tip_label_list[-1], identical, logical(1), tip_label_list[[1]]))
  if(shared){
    clone_info[[paste0(prefix, "_tip_labels_shared")]] <- tip_label_list[[1]]
  }else{
    clone_info[[paste0(prefix, "_tip_labels")]] <- tip_label_list
  }
  clone_info
}

#'\code{dataFrameToValueRows}
#' Like \code{\link{dataFrameToRecords}}, but drops each row's field names,
#' producing a plain JSON array of values per row instead of a JSON object
#' -- for use alongside a column-name vector stored once, separately (the
#' same space-saving idea \code{\link{nexusTreeWithTranslate}} applies to
#' tip names: store the shared "header" once, not repeated per row/tree).
#'
#' Each row is returned as a plain atomic vector when every value in it is
#' numeric and finite: \code{jsonlite::write_json(..., pretty=TRUE)} prints
#' atomic vectors compactly on one line, but a \code{list} of scalars (even
#' all-numeric, unnamed) still gets one element per line. Rows that need
#' string conversion (\code{NaN}/\code{Inf}/\code{-Inf}) -- or, defensively,
#' any row that isn't purely numeric to begin with -- fall back to a
#' \code{list}, so mixed element types are never silently coerced (an atomic
#' vector with any non-numeric value would force every element in that row
#' to character).
#' @param    df  a data.frame, in the fixed column order the row values
#'               should be emitted in
#' @noRd
dataFrameToValueRows <- function(df){
  lapply(seq_len(nrow(df)), function(i){
    row <- unname(as.list(df[i, , drop=FALSE]))
    if(all(vapply(row, is.numeric, logical(1))) &&
      all(vapply(row, is.finite, logical(1)))){
      unlist(row)
    }else{
      lapply(row, function(v) if(is.numeric(v) && !is.finite(v)) as.character(v) else v)
    }
  })
}

#'\code{writeTreesJSON}
#' Experimental. Write trees in AIRR format
#' @param    object   Dowser object from getTrees
#' @param    file     name of .json file
#' @param    repertoire_id  repertoire_id to use
#' @param    check  verify object is equivalent on reloading
#' @param    verbose  print out more info
#' @param    edge_tol tolerance for branch length checks (if check=TRUE)
#' @param    cell     cell id column name in Dowser object
#' @param    heavy    name of heavy chain locus
#' @param    light    names of light chain loci
#' @param   dowser_fields include dowser-specific information? (recommended)
#' @param    nproc    number of cores to use (parallelizes by clone)
#' @details
#' Works with trees built by any of \code{getTrees}'s \code{build} options
#' (\code{"pratchet"}, \code{"pml"}, \code{"igphyml"}, \code{"raxml"}).
#' \code{getTrees(..., build="pml")} trees store the full
#' \code{phangorn::optim.pml} fit in \code{tree$parameters}; since that object
#' isn't JSON-serializable (and isn't meaningful to reconstruct from a file),
#' it's reduced to a flat list of fitted model parameters via
#' \code{pmlParamsToList} before being written out. See that function's
#' documentation for what is kept and why.
#' @export

writeTreesJSON = function(object, file, repertoire_id="sample", check=TRUE, verbose=TRUE,
  edge_tol=1e-8, cell="cell_id", heavy="IGH", light=c("IGK","IGL"), dowser_fields=TRUE,
  nproc=1){

  if(verbose){
    print("Note: internal node numbers/labels not currently preserved.")
  }
  if(inherits(object$trees[[1]], "treedata")){
    timetree <- TRUE
  }else{
    timetree <- FALSE
  }
  clones <- list()
  clones$Clone <- list()
  clones$Rearrangement <- list()
  rearrangement_index <- 1
  for(row in 1:nrow(object)){
    # AIRR standard fields for clone
    clone <- list()
    clone$clone_id <- object$clone_id[row]
    clone$repertoire_id <- repertoire_id
    clone$repertoire_group_id <- NA
    clone$data_processing_id <- NA
    clone$clone_count <- object$seqs[row]
    clone$sequences <- object$data[[row]]@data$sequence_id
    if(dowser_fields){
      # nonstandard but important for dowser
      clone$info <- list()
      if(timetree){
        clone$program_origin <- "dowser-timetree"
      }else{
        clone$program_origin <- "dowser-phylo"
      }
      clone$program_version <- as.character(packageVersion("dowser"))
      #clone$info$locus <- object$data[[row]]@locus
      clone$info$region <- object$data[[row]]@region
      clone$info$numbers <- object$data[[row]]@numbers
      clone$info$phylo_seq <- object$data[[row]]@phylo_seq
      clone$v_call <- object$data[[row]]@v_gene
      clone$j_call <- object$data[[row]]@j_gene
      clone$junction_length <- object$data[[row]]@junc_len
      clone$column_names <- names(object)
    }
    node_class <- "Cell"
    node_index <- "cell_id"
    if(!cell %in% names(object$data[[row]]@data)){
      node_class <- "Rearrangement"
      node_index <- "sequence_id"
      if(object$data[[row]]@phylo_seq == "hlsequence"){
        stop(paste(cell," not found in data object but paired chain data used"))
      }
    }
    clone$clone_class <- node_class
    germline_name <- "Germline"
    germline_qualified <- paste0(germline_name, "-", clone$clone_id)

    if(timetree){
      tr <- object$trees[[row]]
      phy <- tr@phylo
    }else{
      phy <- object$trees[[row]]
    }
    nodes <- length(phy$tip.label) + phy$Nnode
    phy$node.label <- paste0("Node",(length(phy$tip.label)+1):nodes,"-",clone$clone_id)

    phy$tip.label[phy$tip.label == germline_name] = germline_qualified

    airrc <- object$data[[row]]
    germ_seq <- switch(airrc@phylo_seq,
      "sequence" = airrc@germline,
      "lsequence" = airrc@lgermline,
      "hlsequence" = airrc@hlgermline)

    # add info from other columns
    forbidden <- c("region", "numbers", "phylo_seq", "germline",
      "lgermline", "hlgermline","program_origin", "program_version")
    for(n in names(object)){
      if(n %in% c("clone_id","data","locus","seqs","trees")){
        next
      }
      if(n %in% forbidden){
        stop(paste("Forbidden column names:",paste(forbidden, collapse=",")))
      }
      if(dowser_fields){clone$info[[n]] <- object[[n]][[row]]}
    }

    object$data[[row]]@data$tip_order <- 1:nrow(object$data[[row]]@data)
    tipcols <- names(object$data[[row]]@data)
    tipdata <- list()
    for(i in 1:nrow(object$data[[row]]@data)){
      tip <- list()
      for(col in tipcols){
        if(col %in% c("sequence","lsequence","hlsequence")){
          next
        }
        tip[[col]] <- object$data[[row]]@data[i,][[col]]
      }
      tipdata[[tip$sequence_id]] <- tip
    }
    clone$tree <- ape::write.tree(phy)
    if(dowser_fields){
      clone$info$treeinfo <- list()
      for(n in names(phy)){
        if(n %in% c("edge","tip.label","Nnode",
          "edge.length","nodes","node.label", "state")){
          next
        }
        value <- phy[[n]]
        # phangorn::pml fits (getTrees(..., build="pml")) aren't
        # JSON-serializable as-is
        if(n == "parameters" && inherits(value, "pml")){
          value <- pmlParamsToList(value)
        }
        clone$info$treeinfo[[n]] <- value
      }
    }
    ucanode <- ape::getMRCA(phy, tip=phy$tip.label)
    if(sum(phy$tip.label == germline_qualified) > 0){
      germline_node <- which(phy$tip.label == germline_qualified)
    }else{
      germline_node <- NA
      clone$inferred_ancestor <- NA
    }
    # metadata columns for tips
    nodes <- list()
    for(i in 1:(length(phy$tip.label) + phy$Nnode)){
      node <- list()
      if(dowser_fields){
        node$info <- list()
      }
      node$node_class <- node_class
      node$repertoire_id <- repertoire_id
      receptors <- list()
      if(i <= length(phy$tip.label) && (is.na(germline_node) | i != germline_node)){
        node$node_type <- "observed"
        if(dowser_fields){
          node$info$tipdata <- tipdata[[phy$tip.label[i]]]
        }
        seqid <- phy$tip.label[i]
        node$node_id <- seqid
        if(node_class == "Cell"){
            node$cell_id <- object$data[[row]]@data[
            object$data[[row]]@data$sequence_id == seqid,][[cell]]
        }else{
          node$sequence_id <- seqid
        }
      }else{
        node$node_type <- "inferred"
        if(!is.na(germline_node) && i == germline_node){
          node[[node_index]] <- phy$tip.label[i]
        }else{
          node[[node_index]] <- paste0(phy$node.label[i-length(phy$tip.label)])
        }
        node$node_id <- node[[node_index]]
      }
      nodenames <- names(phy$nodes[[i]])
      for(n in nodenames){
        if(n %in% c("id", "sequence")){
          next
        }
        if(dowser_fields){
          node$info[[n]] <- unlist(phy$nodes[[i]][n])
        }
      }
      node$number <- i
      node$parent <- phy$edge[phy$edge[,2] == i,1]
      node$children <- phy$edge[phy$edge[,1] == i,2]

      # need to update this by getting getNodeSeq to work on time trees
      # add name field to treedata@info
      # return null if no nodes list unless requesting tip sequence
      seqa <- getNodeSeq(object, node=i, tree=phy, gaps=TRUE)
      seq <- getNodeSeq(object, node=i, tree=phy, gaps=FALSE)
      if(!is.na(germline_node)){
        ga <- getNodeSeq(object, node=germline_node, tree=phy, gaps=TRUE)
      }else{
        ga <- NA
      }
      for(loci_i in 1:length(seq)){
        receptor <- list()
        loci <- names(seq)[loci_i]
        addon <- ""
        if(node_class == "Cell"){
          receptor$cell_id <- node$cell_id
          addon <- paste0("-",loci)
        }
        receptor$sequence_id <- paste0(node[[node_index]], addon)
        receptor$sequence_alignment <- seqa[loci_i]
        receptor$sequence <- seq[loci_i]
        receptor$germline_alignment <- ga[loci_i]
        receptor$locus <- loci
        receptor$rev_comp <- NA
        receptor$productive <- NA
        receptor$v_call <- NA
        receptor$d_call <- NA
        receptor$j_call <- NA
        receptor$junction <- NA
        receptor$junction_aa <- NA
        receptor$junction_length <- NA
        receptor$v_cigar <- NA
        receptor$d_cigar <- NA
        receptor$j_cigar <- NA
        receptor$rearrangement_type <- node$node_type
        clones$Rearrangement[[rearrangement_index]] <- receptor
        rearrangement_index <- rearrangement_index + 1
      }
      if(!is.null(phy$state)){
        node$state_vector_value <- phy$state[i]
      }
      #node$node_id <- node[[node_index]]
      nodes[[i]] <- node

      if(!is.na(germline_node) && node$number == germline_node){
        clone$inferred_ancestor <- node$node_id
      }
    }
    clone$nodes <- nodes

    if(timetree){
      clone$info$treetext <- tr@treetext
      main_tip_labels <- phy$tip.label
      if(!is.na(germline_node)){
        main_tip_labels[main_tip_labels == germline_name] <- germline_qualified
      }
      clone$info$tip_labels <- main_tip_labels
      #clone$info$root_edge <- tr@phylo$root.edge
      clone$info$beast_parameters <- dataFrameToRecords(tr@info$parameters)
      # Full posterior (full_posterior=TRUE), when present -- see this
      # function's @details. tr@info$trees_posterior (the raw sampled trees)
      # and tr@info$trees_with_traits_posterior (a separate, typically
      # thinned, set of posterior trees carrying ancestral trait
      # reconstruction in @data) are two independent, optional fields --
      # any combination of them (and parameters_posterior) may be present.
      # Each is stored the same way, via writePosteriorTreeList().
      for(prefix in c("trees_posterior", "trees_with_traits_posterior")){
        posterior_trees <- tr@info[[prefix]]
        if(!is.null(posterior_trees)){
          clone$info <- writePosteriorTreeList(clone$info, posterior_trees,
            prefix, germline_name, germline_qualified, !is.na(germline_node))
        }
      }
      if(!is.null(tr@info$parameters_posterior)){
        #wide_pp <- tidyr::pivot_wider(tr@info$parameters_posterior,
        #  names_from="parameter", values_from="value")
        wide_pp <- tr@info$parameters_posterior
        clone$info$parameters_posterior_columns <- names(wide_pp)
        clone$info$parameters_posterior <- dataFrameToValueRows(wide_pp)
      }
    }

    clones$Clone[[row]] <- clone

  }
  if(grepl(".json$",file)){
    # digits=NA keeps small floating points from getting truncated
    jsonlite::write_json(clones, file, pretty=TRUE, simplifyVector=TRUE,
      auto_unbox = TRUE, digits=NA)
  }else{
    stop("file extension not recognized, must end in .json")
  }
  # verify object was stored faithfully
  if(check){
    if(verbose){
      print("Loading object to check consistency")
    }
    nobject <- readTreesJSON(file, heavy=heavy, light=light, verbose=verbose, nproc=nproc)
    # dowserObjectEquivalent understands pml fits directly (reducing them via
    # pmlParamsToList itself), so no pre-normalization is needed here.
    validate <- dowserObjectEquivalent(object, nobject, verbose=verbose, 
      edge_tol=edge_tol, dowser_fields=dowser_fields, nproc=nproc)
  }
}




#' \code{readTreesJSON}
#' Experimental. Read trees from JSON/AIRR format from Dowser
#' @param    file   .json file
#' @param    heavy    name of heavy chain locus
#' @param    light    names of light chain loci
#' @param    verbose  how much info to print
#' @param    edge_tol tolerance for branch length checks (if check=TRUE)
#' @param    nproc    number of cores to use (parallelizes by clone)
#' @details
#' Reads files written by \code{\link{writeTreesJSON}}, including trees built
#' with any \code{getTrees} \code{build} option.
#' @export
readTreesJSON = function(file, heavy="IGH", light=c("IGK","IGL"), 
  verbose=TRUE, edge_tol=1e-8, nproc=1){

  rclones <- jsonlite::read_json(file)
  program_origin <- rclones$Clone[[1]]$program_origin

  output <- dplyr::tibble()
  outtrees <- list()
  outdata <- list()
  rearrangements <- rclones$Rearrangement

  if("program_origin" %in% names(rclones$Clone[[1]]) &&
    rclones$Clone[[1]]$program_origin == "dowser-phylo"){
    dowser <- TRUE
    timetree <- FALSE
  }else if("program_origin" %in% names(rclones$Clone[[1]]) &&
    rclones$Clone[[1]]$program_origin == "dowser-timetree"){
    dowser <- TRUE
    timetree <- TRUE
  }else{
    if(verbose)warning("Some features limited by using non-Dowser JSON")
    dowser <- FALSE
    timetree <- FALSE
  }

  node_class <- rclones$Clone[[1]]$clone_class
  if(node_class == "Cell"){
    node_index <- "cell_id"
  }else{
    node_index <- "sequence_id"
  }
  rearrangement_names <- sapply(rearrangements, function(x)x[[node_index]])
  for(ci in 1:length(rclones$Clone)){
    clone <- rclones$Clone[[ci]]
    seqs <- clone$clone_count[[1]]
    rphy <- ape::read.tree(text=clone$tree)
    nodes <- clone$nodes

    loci <- c()
    numbers <- c()
    column_names <- c()
    ancestor_node_id <- ifelse(is.null(clone$inferred_ancestor),
      NA,clone$inferred_ancestor)

    rnodes <- list()
    tips <- list()
    tipcount <- 1
    for(i in 1:length(nodes)){
      node <- nodes[[i]]
      #node_num =  node$info$number[[1]]
      if(node$node_id %in% rphy$tip.label){
        new_num <- which(rphy$tip.label == node$node_id)
      }else if(node$node_id %in% rphy$node.label){
        new_num <- which(rphy$node.label == node$node_id) + length(rphy$tip.label)
      }else{
        stop(paste(node$node_id, "not found in tree"))
      }
      new_node <- list()
      for(name in names(node$info)){
        if(name %in% c("number","parent","children","tipdata")){
          next
        }else if(name == "state_vector_value"){
          if(is.null(rphy$state)){
            rphy$state <- rep(NA, length=length(nodes))
          }
          rphy$state[new_num] <- unlist(node$info[[name]])
        }else{
          new_node[[name]] <- node$info[[name]]
        }
      }
      receptors <- rearrangements[rearrangement_names == node[[node_index]]]
      receptor_loci <- sapply(receptors, function(x)x$locus)
      loci_order <- order(receptor_loci)
      if(node_class == "Cell" & receptor_loci[loci_order[1]] != heavy){
        stop(paste(heavy, " not in first ordered locus position"))
      }

      # assumes reeptors are in same order as locus object
      seq <- ""
      for(lo in loci_order){
        seq <- paste0(seq, receptors[[lo]]$sequence)
        if(i == 1){
          loci <- c(loci, rep(receptors[[lo]]$locus,
            nchar(receptors[[lo]]$sequence)))
          if(!dowser){
            numbers <- c(numbers, 1:nchar(receptors[[lo]]$sequence))
          }
        }
      }
      if(i == 1){
        unique_loci <- unique(loci)
        if(dowser && !is.null(clone$info$phylo_seq)){
          # Dowser-written files record the true phylo_seq; use it
          # directly instead of re-deriving it from the locus tag
          phylo_seq <- clone$info$phylo_seq[[1]]
          column_names <- unlist(clone$column_names)
        }else if(heavy %in% loci && length(intersect(light,unique_loci)) > 0){
          phylo_seq <- "hlsequence"
        }else if(unique_loci[1] == heavy && length(unique_loci) == 1){
          phylo_seq <- "sequence"
        }else if(length(intersect(light,unique_loci)) > 0){
          phylo_seq <- "lsequence"
        }else{
          warning(paste("Could assign phylo_seq based on locus column, using 'sequence'",
            paste(unique_loci, collapse=",")))
          phylo_seq <- "sequence"
        }
      }
      if(nchar(seq) != length(loci)){
        if(!timetree)stop("seq and loci vector not the same length!")
      }
    if(seq == ""){seq <- NA}
      new_node$sequence <- seq
      node$sequence <- seq
      rnodes[[new_num]] <- new_node
      if(node$node_id %in% rphy$tip.label){
        tips[[tipcount]] <- node
        tipcount <- tipcount + 1
      }
    }
    rphy$nodes <- rnodes
    rphy$node.label <- NULL
    rphy$name <- clone$clone_id

    for(n in names(clone$info$treeinfo)){
      if(n == "state"){
        next
      }
      ni <- clone$info$treeinfo[[n]]
      if(length(ni) > 1){
        ni <- lapply(ni, function(x)unlist(x))
        if(n == "parameters" && !is.null(clone$info$treeinfo$tree_method) &&
          grepl("optim\\.pml", clone$info$treeinfo$tree_method)){
          tip_seqs <- stats::setNames(
            sapply(tips, function(tp) tp$sequence),
            sapply(tips, function(tp) tp$node_id))
          tree_for_pml <- rphy
          #only mess with the germline name if coming from Dowser
          if(dowser && !is.na(ancestor_node_id)){
            names(tip_seqs)[names(tip_seqs) == ancestor_node_id] <- "Germline"
            tree_for_pml$tip.label[tree_for_pml$tip.label == ancestor_node_id] <- "Germline"
          }
          ni <- listToPmlParams(ni, tree=tree_for_pml, seqs=tip_seqs)
        }
        rphy[[n]] <- ni
      }else{
        rphy[[n]] <- ni[[1]]
      }
    }

    #extract data from the tips
    data <- list()
    germline <- NA
    for(i in 1:length(tips)){
      tip <- tips[[i]]
      if(dowser){
        row <- tip$info$tipdata
      }else{
        row <- list()
        row$sequence_id <- tip$node_id
      }
      row[[phylo_seq]]<- tip$sequence

      if(is.na(ancestor_node_id) || tip$node_id != ancestor_node_id){
        data[[i]] <- row
      }else{
        germline <- row[[phylo_seq]]
      }
    }
    bdata <- dplyr::bind_rows(data)
    alignment_width <- unique(nchar(bdata[[phylo_seq]]))
    if(length(alignment_width) > 1){
      stop(paste(rclones$Clone[[ci]]$clone_id, "data sequences different length"))
    }
    if(is.na(germline)){
      germline <- paste(rep("N",times=alignment_width), collapse="")
    }

    if(dowser){
      bdata <- bdata[order(bdata$tip_order),]
      bdata <- dplyr::select(bdata, -"tip_order")
      if(!is.na(ancestor_node_id)){
        rphy$tip.label[rphy$tip.label == ancestor_node_id] <- "Germline"
      }
    }

    lgermline <- ""
    hlgermline <- ""
    if(!dowser){
      region <- rep("N",times=alignment_width)
    }
    if(phylo_seq == "hlsequence"){
      hlgermline <- germline
      lgermline <- ""
      germline <- ""
      if(!dowser){
        region <- rep("N",times=alignment_width)
      }
    }else if(phylo_seq == "lsequence"){
      lgermline <- germline
      germline <- ""
      hlgermline <- ""
      if(!dowser){
        region <- rep("N",times=alignment_width)
      }
    }

    if(dowser){
      region <- unlist(clone$info$region)
      numbers <- unlist(clone$info$numbers)
    }

    outclone <- new("airrClone",
          data=bdata,
          clone=as.character(clone$clone_id[[1]]),
          germline=germline,
          lgermline=lgermline,
          hlgermline=hlgermline,
          v_gene=ifelse("v_call" %in% names(clone),unlist(clone$v_call),""),
          j_gene=ifelse("j_call" %in% names(clone),unlist(clone$j_call),""),
          junc_len=ifelse("junction_length" %in% names(clone),clone$junction_length[[1]],0),
          locus=loci,
          region=region,
          numbers=numbers,
          phylo_seq=phylo_seq)

    outdata[[ci]] <- outclone
    
    if(!timetree){
      outtrees[[ci]] <- rphy
    }else{
      # specialized functions for reading in treedata objects
      tr <- treeio::read.beast.newick(textConnection(clone$info$treetext[[1]]))
      tr@phylo$tip.label <- unname(sapply(tr@phylo$tip.label, 
        function(x)clone$info$tip_label[[as.numeric(x)]]))
      if(!is.na(ancestor_node_id)){
        tr@phylo$tip.label[tr@phylo$tip.label == ancestor_node_id] <- "Germline"
      }

      # check and copy over attributes that were read in from the nodes list
      # this might seem redundant but ensures that two versions of the same
      # tree are consistent.
      # node_map, rphy node = node_map[tr@phylo node]
      if(!treesEquivalent(rphy, tr@phylo, edge_tol=edge_tol)){
        stop(paste(clone$clone_id, "tree and treetext trees not the same"))
      }
      node_map <- mapSubtrees(tr@phylo, rphy)
      tr@phylo$nodes <- lapply(1:length(rphy$nodes), function(x)rphy$nodes[[node_map[x]]])
      for(n in names(rphy)){
        if(n %in% c("edge", "edge.length", "Nnode", "tip.label", "nodes")){
          next
        }else{
          tr@phylo[[n]] <- rphy[[n]]
        }
      }
      #tr@phylo$root.edge <- if(is.null(clone$info$root_edge)) 0 else clone$info$root_edge[[1]]

      # read in beast parameter estimates
      beast_params <- clone$info$beast_parameters
      if(!is.null(beast_params)){
        params_df <- dplyr::bind_rows(lapply(beast_params, function(row){
          as.data.frame(lapply(row, function(x) as.character(x[[1]])), stringsAsFactors=FALSE)
        }))
        # coerce every column but "item" back to numeric
        for(col in names(params_df)){
          if(col != "item"){
            params_df[[col]] <- suppressWarnings(as.numeric(params_df[[col]]))
          }
        }
        tr@info$parameters <- params_df
      }

      # read in tree posterior fields
      for(prefix in c("trees_posterior", "trees_with_traits_posterior")){
        if(prefix %in% names(clone$info)){
          tips_shared <- clone$info[[paste0(prefix, "_tip_labels_shared")]]
          tips_raw <- clone$info[[paste0(prefix, "_tip_labels")]]
          raw <- clone$info[[prefix]]
          tr@info[[prefix]] <- parallel::mclapply(seq_along(raw), function(k){
            tip_labels_k <- if(!is.null(tips_shared)){
              unlist(tips_shared)
            }else{
              unlist(tips_raw[[k]])
            }
            pt <- treeio::read.beast.newick(textConnection(raw[[k]]))
            pt@phylo$tip.label <- unname(sapply(pt@phylo$tip.label, 
              function(x)tip_labels_k[[as.numeric(x)]]))
            if(!is.na(ancestor_node_id)){
              pt@phylo$tip.label[pt@phylo$tip.label == ancestor_node_id] <- "Germline"
            }
            pt
          }, mc.cores=nproc)
        }
      }

      # read in parameters posterior
      # [[ ]], not $: consistent with the posterior-tree fields above.
      params_post_cols <- clone$info[["parameters_posterior_columns"]]
      params_post_raw <- clone$info[["parameters_posterior"]]
      if(!is.null(params_post_raw)){
        # Stored wide, as column names (once) + one unnamed value array per
        # Sample -- reattach the column names to each row first (same
        # all-character-first fix as beast_params above; every column is
        # numeric here, no "item"/"parameter" name column to exclude, since
        # parameter names are now the column headers), then pivot back to the
        # long Sample/parameter/value shape tr@info$parameters_posterior
        # actually has (matching tidyr::gather()'s output in
        # R/TimeTreesFunctions.R), so the reconstructed object is a faithful
        # round-trip, not just equivalent data in a different shape.
        cols <- unlist(params_post_cols)
        params_post_wide <- dplyr::bind_rows(lapply(params_post_raw, function(row){
          vals <- vapply(row, as.character, character(1))
          as.data.frame(as.list(stats::setNames(vals, cols)), stringsAsFactors=FALSE)
        }))
        for(col in names(params_post_wide)){
          params_post_wide[[col]] <- suppressWarnings(as.numeric(params_post_wide[[col]]))
        }
        tr@info$parameters_posterior <- params_post_wide
        #TODO: change readBEAST and functions to keep these as wide-format to save space
        #tr@info$parameters_posterior <- tidyr::pivot_longer(params_post_wide,
        #  cols=-"Sample", names_to="parameter", values_to="value")
      }
      outtrees[[ci]] <- tr
    }

    # store info in other columns
    temp <- dplyr::tibble(clone_id=clone$clone_id[[1]], seqs=nrow(outclone@data),
      locus=paste0(sort(unique(loci)),collapse=","))
    for(n in names(clone$info)){
      if(!n %in% c("region", "numbers", "phylo_seq", "germline", "lgermline",
        "hlgermline", "trees", "data", "clone_id", "seqs", "locus", "program_origin", "treeinfo",
        "program_version", "column_names", 
        "treetext", "tip_labels", "root_edge", "beast_parameters",
        "trees_posterior", "trees_posterior_tip_labels", "trees_posterior_tip_labels_shared",
        "trees_with_traits_posterior", "trees_with_traits_posterior_tip_labels",
        "trees_with_traits_posterior_tip_labels_shared",
        "parameters_posterior", "parameters_posterior_columns")){
        #print(n)
        ni <- clone$info[[n]]
        if(timetree && n %in% c("parameters","skyline")){
          ni <- as.data.frame(dplyr::bind_rows(ni))
          temp[[n]] <- list(ni)
        }else if(length(ni) > 1){
          ni <- lapply(ni, function(x)unlist(x))
          temp[[n]] <- list(ni)
        }else{
          temp[[n]] <- ni[[1]]
        }
      }
    }
    output <- dplyr::bind_rows(output, temp)
  }
  output$data <- outdata
  output$trees <- outtrees
  onames <- names(output)

  if(dowser){
    if(length(setdiff(onames, column_names)) > 0){
        print(paste(onames))
        print(paste(column_names))
        stop("mismatch in object column names")
      }
      output <- output[,column_names]
  }else{
    standard <- c("clone_id", "data", "seqs", "locus", "trees")
    nonstandard <- setdiff(onames, standard)
    output <- dplyr::select(output, dplyr::any_of(standard), dplyr::any_of(nonstandard))
  }
  output
}

#'\code{pmlParamsEqual}
#' Compare two sets of \code{phangorn::pml} fit parameters for equality,
#' within a numeric tolerance.
#'
#' Accepts either a full \code{pml} fit object (as stored in \code{tree$parameters}
#' after \code{getTrees(..., build="pml")}) or the flat list produced by
#' \code{\link{pmlParamsToList}} (as stored after a JSON round-trip) on either
#' side. Also works with RAxML and IgPhyML parameter lists
#'
#' When \emph{both} sides are full \code{pml} objects with a real
#' \code{$data}/\code{$tree} (i.e. both went through
#' \code{\link{listToPmlParams}}'s reconstruction -- \code{igphyml}/
#' \code{raxml} parameter lists never have these), the real alignment and
#' real tree are also compared: alignment content via \code{as.character()}
#' on the \code{phyDat} (tip order-independent, matched by name), and the
#' tree via pairwise tip-to-tip (patristic) distances from
#' \code{ape::cophenetic.phylo()} rather than a direct edge-by-edge
#' comparison -- topologically identical trees can have their internal nodes
#' numbered/ordered differently (e.g. after an \code{ape::unroot()}), and
#' cophenetic distances are invariant to that while still fully capturing
#' both topology and branch lengths.
#' @param    pa        first \code{parameters} value (\code{pml} object or plain list)
#' @param    pb        second \code{parameters} value (\code{pml} object or plain list)
#' @param    tol       relative tolerance (passed to \code{all.equal}) for comparing
#'                     fitted scalar/vector values
#' @param    edge_tol  absolute tolerance for comparing the real tree's pairwise tip
#'                     distances, when both sides have a tree
#' @noRd
pmlParamsEqual <- function(pa, pb, tol=1e-3, edge_tol=1e-8){
  pa_full <- pa
  pb_full <- pb
  if(inherits(pa, "pml")){ pa <- pmlParamsToList(pa) }
  if(inherits(pb, "pml")){ pb <- pmlParamsToList(pb) }
  if(!setequal(names(pa), names(pb))){
    return(FALSE)
  }
  for(nm in names(pa)){
    # unname(): a named vector (e.g. pmlParamsToList's labeled bf/Q) survives
    # a JSON round-trip as a plain unnamed array, so names must be ignored
    # here or a real value match still reports as unequal.
    #va <- uname(unlist(pa[[nm]]))
    #vb <- uname(unlist(pb[[nm]]))
    va <- unlist(pa[[nm]])
    vb <- unlist(pb[[nm]])
    if(length(va) != length(vb)){
      return(FALSE)
    }
    if(is.numeric(va) && is.numeric(vb)){
      if(!isTRUE(all.equal(va, vb, tolerance=tol))){
        return(FALSE)
      }
    }else if(any(as.character(va) != as.character(vb))){
      return(FALSE)
    }
  }

  # Real alignment, when both sides have one (only true full pml fits do --
  # igphyml/raxml parameter lists never carry data/tree). 
  # [[ ]], not $ avoid partial tree and tree_length match with raxml/igphyml
  data_a <- pa_full[["data"]]
  data_b <- pb_full[["data"]]
  if(is.null(data_a) != is.null(data_b)){
    return(FALSE)
  }
  if(!is.null(data_a)){
    da <- as.character(data_a)
    db <- as.character(data_b)
    if(!setequal(rownames(da), rownames(db))){
      return(FALSE)
    }
    db <- db[rownames(da), , drop=FALSE]
    if(!identical(dim(da), dim(db)) || sum(da != db) != 0){
      return(FALSE)
    }
  }

  tree_a <- pa_full[["tree"]]
  tree_b <- pb_full[["tree"]]
  if(is.null(tree_a) != is.null(tree_b)){
    return(FALSE)
  }
  if(!is.null(tree_a)){
    if(!setequal(tree_a$tip.label, tree_b$tip.label)){
      return(FALSE)
    }
    ca <- ape::cophenetic.phylo(tree_a)
    cb <- ape::cophenetic.phylo(tree_b)[rownames(ca), rownames(ca)]
    if(!isTRUE(all.equal(ca, cb, tolerance=edge_tol, check.attributes=FALSE))){
      return(FALSE)
    }
  }

  return(TRUE)
}

#' Get a map of nodes from one 
#' @param    treea  First phylo object
#' @param    treeb  Second phylo object
#' @details 
#' Returns a vector in which each position is the numerical
#' node number in treea, and the value is the corresponding node number in treeb
#' So node 1 in treea corresponds to result[1] in treeb
#' @export
mapSubtrees = function(treea, treeb){
  nodesb <- c() #map of nodes from a to b
  tipsa <- treea$tip.label
  tipsb <- treeb$tip.label
  if(!setequal(tipsa, tipsb)){
    diffs <- c(setdiff(tipsb, tipsa), setdiff(tipsa, tipsb))
    fdiffs <- diffs[!grepl("Germline", diffs)]
    if(length(fdiffs) == 0){
      treea$tip.label[treea$tip.label %in% diffs] = "Germline"
      treeb$tip.label[treeb$tip.label %in% diffs] = "Germline"
    }else{
      stop(paste("tips not the same"))
    }
  }
  all_subtrees_a <- lapply(1:(length(treea$tip.label) + treea$Nnode), function(x)getSubTaxa(x, treea))
  all_subtrees_b <- lapply(1:(length(treeb$tip.label) + treeb$Nnode), function(x)getSubTaxa(x, treeb))
  if(length(all_subtrees_a) != length(all_subtrees_b)){
    stop(paste("trees don't have the same number of subtrees"))
  }
  # check if trees have the same subtrees with corresponding edge lengths and sequences
  for(sta in 1:length(all_subtrees_a)){
    sa <- all_subtrees_a[[sta]]
    match <- -1
    for(stb in 1:length(all_subtrees_b)){
      if(setequal(sa, all_subtrees_b[[stb]])){
        match <- stb
      }
    }
    if(match < 0){
      stop(paste("subtrees don't match"))
      treecheck = FALSE
    }
    nodesb[sta] <- match
  }
  nodesb
}


#'\code{dataColumnEqual}
#' Compare one \code{@data} column between two BEAST-annotated trees,
#' tolerant of a floating-point string-formatting artifact: BEAST numeric
#' annotations are stored as \emph{character} strings (not \code{numeric}),
#' and \code{treeio::read.beast()} is not idempotent on them -- re-parsing a
#' tree's own \code{@treetext} a second time can shift the last
#' significant digit of a numeric annotation's string representation (e.g.
#' \code{"0.51802396741998"} vs \code{"0.518023967419979"} for the exact
#' same underlying value), confirmed by re-parsing a real object's own
#' \code{@treetext} independently of any write/read round trip. Columns
#' that parse fully as numeric are compared with tolerance; anything else
#' (e.g. categorical trait columns like \code{location}) is compared
#' exactly, as before.
#' @param    va, vb     the two columns to compare
#' @param    tolerance  numeric tolerance for columns that parse as numeric
#' @noRd
dataColumnEqual <- function(va, vb){
  if(is.character(va) && is.character(vb)){
    na_num <- suppressWarnings(as.numeric(va))
    nb_num <- suppressWarnings(as.numeric(vb))
    if(!anyNA(na_num) && !anyNA(nb_num)){
      return(isTRUE(all.equal(na_num, nb_num, check.attributes=FALSE)))
    }else{
      # KBH
      # if same indexes are NA, compare values that aren't
      if(sum(is.na(na_num) != is.na(nb_num)) == 0){
        na_index <- !is.na(na_num)
        return(isTRUE(all.equal(na_num[na_index], nb_num[na_index], 
          tolerance=tolerance, check.attributes=FALSE)))
      }else{
        return(FALSE)
      }
    }
  }
  isTRUE(all.equal(va, vb))
}

#' Check whether two tree objects are equivalent
#' @param    obja  First phylo or treedata object
#' @param    objb  Second phylo or treedata object
#' @param    edge_tol tolerance for branch length checks (if check=TRUE)
#' @param    numbering_match require internal node numbers to match?
#' @param    clonesa  Dowser clones object associated with obja (if check_extended=TRUE)
#' @param    clonesb  Dowser clones object associated with objb (if check_extended=TRUE)
#' @param    check_extended Also check node sequences and state vector? Requires clonesa/b to be supplied
#' @param    gaps  check sequences wth IMGT gaps with check_extended?
#' @details For treedata objects, check both @phylo and @data
#' @export
treesEquivalent = function(obja, objb, edge_tol=1e-8, numbering_match=FALSE,
  clonesa=NULL, clonesb=NULL, check_extended=FALSE, gaps=TRUE){
  
  treecheck <- TRUE
  if(check_extended){
    if(is.null(clonesa) || is.null(clonesb)){
      stop("clonesa and clonesb must be provided if check_extended=TRUE")
    }
  }
  a_is_timetree <- inherits(obja, "treedata")
  b_is_timetree <- inherits(objb, "treedata")
  if(a_is_timetree != b_is_timetree){
    stop("Cannot compare a time-tree object against a non-time-tree object")
  }
  if(a_is_timetree){
    timetree <- TRUE
    treea <- obja@phylo
    treeb <- objb@phylo
  }else{
    timetree <- FALSE
    treea <- obja
    treeb <- objb
  }

  node_map <- tryCatch(mapSubtrees(treea, treeb), error=function(e)e)
  if(inherits(node_map, "error")){
    warning(node_map)
    return(FALSE)
  }
  # check if trees have the same subtrees with corresponding edge lengths and sequences
  for(nodea in 1:length(node_map)){
    nodeb <- node_map[nodea]
    ea <- round(treea$edge.length[treea$edge[,2] == nodea], digits=11)
    eb <- round(treeb$edge.length[treeb$edge[,2] == nodeb], digits=11)
    if(length(ea) > 0 || length(eb) > 0){
      if(abs(ea - eb) > edge_tol){
        warning(paste( "Edges not within edge_tol", ea, eb, nodea, nodeb))
        treecheck = FALSE
      }
    }
    if(check_extended){
      seqa <- getNodeSeq(clonesa, node=nodea, tree=treea, gaps=gaps)
      seqb <- getNodeSeq(clonesb, node=nodeb, tree=treeb, gaps=gaps)
      if(!is.null(seqa) && !is.null(seqb)){
        treecheck = isTRUE(all.equal(seqa, seqb))
        if(!treecheck){
          warning(paste("sequences not identical", nodea, nodeb))
        }
      }
      if(!is.null(treea$state) || !is.null(treeb$state)){
        if(treea$state[nodea] != treeb$state[nodeb]){
          warning(paste("node states not equal",treea$state[nodea],
            treea$state[nodeb], nodea, nodeb))
        }
      }
    }
  }

  if(timetree){
     da <- obja@data
     db <- objb@data
     # node numbers aren't necessarily the same, so map one to the other
     da$nodeb <- as.character(node_map[as.numeric(da$node)])
     da <- da[order(da$nodeb),]
     db <- db[order(db$node),]
     if(!identical(da$node, db$node) && numbering_match){
       warning(paste("@data node numbering not the same"))
       treecheck <- FALSE
     }
     #da <- select(da, -!!rlang::sym("nodeb"))
     da <- da[,names(da) != "nodeb"]
     if(!setequal(names(da), names(db))){
       warning(paste("@data columns not the same"))
       treecheck <- FALSE
     }
     for(col in names(da)){
       if(col == "node"){
        next
       }
       if(!dataColumnEqual(da[[col]], db[[col]])){
         warning(paste("@data column", col, "not the same"))
         treecheck <- FALSE
       }
     }
  }

  treecheck
}



#'\code{dowserObjectEquivalent}
#' Experimental. Check if two Dowser objects are equivalent
#' @param    obj1  First Dowser object
#' @param    obj2  Second Dowser object
#' @param    verbose  print out more info
#' @param    edge_tol tolerance for branch length checks (if check=TRUE)
#' @param   dowser_fields check dowser-specific fields and gapped sequences?
#' @param    nproc number of cores to use
#' @details
#' In addition to the existing tree topology, edge length, sequence, and
#' data slot checks, this also verifies \code{tree$parameters} when present
#' -- including \code{build="pml"} trees.
#' @export
dowserObjectEquivalent = function(obj1, obj2, verbose=TRUE, edge_tol=1e-8,
  dowser_fields=TRUE, nproc=1){

  a_is_timetree <- inherits(obj1$trees[[1]], "treedata")
  b_is_timetree <- inherits(obj2$trees[[1]], "treedata")
  if(a_is_timetree != b_is_timetree){
    stop("Cannot compare a time-tree object against a non-time-tree object")
  }
  timetree <- FALSE
  if(a_is_timetree){
    timetree <- TRUE
  }

  a <- obj1
  b <- obj2
  a <- a[order(a$clone_id),]
  b <- b[order(b$clone_id),]

  if(nrow(a) != nrow(b)){
    stop("unequal row numbers")
  }
  if(sum(a$clone_id != b$clone_id)){
    stop("different clone ids")
  }
  if(sum(a$seqs != b$seqs) != 0){
    stop("different seq numbers")
  }
  if(sum(a$locus != b$locus) != 0){
    stop("different locus columns")
  }
  if(dowser_fields && sum(names(a) != names(b)) > 0){
    stop("different column names")
  }

  #for(r in 1:nrow(a)){
  checks <- lapply(1:nrow(a),function(r)tryCatch({
    treecheck <- 0

    # check non-tree columns
    if(dowser_fields){
      for(n in names(a)){
        if(n %in% c("data", "trees")){
          next
        }
        vala <- a[[n]][[r]]
        valb <- b[[n]][[r]]
        if(inherits(vala,"data.frame")){
          check <- isTRUE(all.equal(vala, valb))
          if(!check){
            stop(paste(n, r, "columns not equal"))
          }
        }
      }
    }

    # check trees
    if(timetree){
      treea <- a$trees[[r]]@phylo
      treeb <- b$trees[[r]]@phylo
      treecheckfunc <- treesEquivalent(a$trees[[r]], b$trees[[r]], edge_tol, 
        clonesa=obj1, clonesb=obj2, check_extended=TRUE, gaps=dowser_fields)
    }else{
      treea <- a$trees[[r]]
      treeb <- b$trees[[r]]
      treecheckfunc <- treesEquivalent(treea, treeb, edge_tol, 
        clonesa=obj1, clonesb=obj2, check_extended=TRUE, gaps=dowser_fields)
    }
    if(!treecheckfunc){
      stop(paste(r, "trees not equivalent (see warnings)"))
    }else{
      treecheck <- treecheck + 1
    }
    # check remaining tree info
    if(dowser_fields && !setequal(names(treea), names(treeb))){
      namediff <- setdiff(names(treea), names(treeb))
      if(length(namediff) > 1 || namediff != "node.label"){
        print(namediff)
        stop(paste(a$clone_id[r], "Tree names not equal"))
      }
      if("node.label" %in% c(names(treea),names(treeb))){
        if(r == 0){
          warning("node.label not currently checked, or preserved by writeTreesJSON")
        }
      }
    }else{
      treecheck <- treecheck + 1
    }
    for(n in names(treea)){
      if(n %in% c("edge","tip.label","Nnode","edge.length","nodes","node.label")){
        next
      }
      if(n == "parameters" && dowser_fields){
        if(is.null(treea[[n]]) && is.null(treeb[[n]])){
          treecheck <- treecheck + 1
        }else if(is.null(treea[[n]]) || is.null(treeb[[n]])){
          stop(paste(a$clone_id[r], n, "not the same"))
        }else if(!pmlParamsEqual(treea[[n]], treeb[[n]], edge_tol=edge_tol)){
          stop(paste(a$clone_id[r], n, "not the same"))
        }else{
          treecheck <- treecheck + 1
        }
        next
      }
      null <- FALSE
      if(sum(!is.na(treea[[n]])) || sum(!is.null(treea[[n]]))){
        if(sum(!is.na(treeb[[n]])) || sum(!is.null(treeb[[n]]))){
          null <- TRUE
        }
      }
      if(!null && sum(treea[[n]] != treeb[[n]]) != 0){
        stop(paste(a$clone_id[r], n, "not the same"))
      }else{
        treecheck <- treecheck + 1
      }
    }

    # check data
    da <- a$data[[r]]
    db <- b$data[[r]]
    da@data <- da@data[order(da@data$sequence_id),]
    db@data <- db@data[order(db@data$sequence_id),]
    for(n in slotNames(da)){
      if(n == "data"){
        dataa = slot(da,n)
        datab = slot(db,n)

        diffnames = setdiff(names(dataa), names(datab))
        if(dowser_fields && sum(!diffnames %in% c("sequence","hlsequence",
          "lsequence","tip_order"))){
            stop(paste("data names don't match",paste(diffnames, collapse=",")))
        }

        namechecks <- union(names(dataa), names(datab))
        if(!dowser_fields){
          namechecks <- intersect(names(dataa), names(datab))
          if(!"sequence_id" %in% namechecks || !da@phylo_seq %in% namechecks){
            stop("Minimal data columns not present in both")
          }
        }else{
          namechecks <- namechecks[!namechecks %in% c("sequence","hlsequence",
            "lsequence","tip_order","sequence_id")]
          namechecks <- c(namechecks, "sequence_id", da@phylo_seq)
        }
        for(name in namechecks){
          #if(sum(dataa[[name]] != datab[[name]]) != 0){
          if(!isTRUE(all.equal(dataa[[name]], datab[[name]]))){
            stop(paste(a$clone_id[r], "data", name, "not the same"))
          }
        }
      }else if(n == "germline" & da@phylo_seq != "sequence"){
        next
      }else if(n == "lgermline" & da@phylo_seq != "lsequence"){
        next
      }else if(n == "hlgermline" & da@phylo_seq != "hlsequence"){
        next
      }else if(!dowser_fields && n %in% c("v_gene","j_gene","junc_len",
        "numbers", "region")){
        next
      }else if(sum(slot(da,n) != slot(db, n)) != 0){
        stop(paste(a$clone_id[r], n, "slot is not equal"))
      }else{
        treecheck <- treecheck + length(slot(da,n))
      }
    }

    #time tree object checks
    if(timetree){
     dtreea <- a$trees[[r]]
     dtreeb <- b$trees[[r]]
     
     # BEAST parameter log 
     pa <- dtreea@info$parameters
     pb <- dtreeb@info$parameters
     if(is.null(pa) != is.null(pb)){
       stop(paste(a$clone_id[r], "beast parameters presence not the same"))
     }
     if(!is.null(pa)){
       if(!setequal(pa$item, pb$item)){
         stop(paste(a$clone_id[r], "beast parameter items not the same"))
       }
       pb <- pb[match(pa$item, pb$item),]
       for(col in setdiff(names(pa), "item")){
         va <- pa[[col]]; vb <- pb[[col]]
         ok <- isTRUE(all.equal(va, vb)) ||
           all(is.na(va) == is.na(vb) & (is.na(va) | va == vb))
         if(!ok){
           stop(paste(a$clone_id[r], "beast parameter", col, "not the same"))
         }
       }
     }
     treecheck <- treecheck + 1
  
      # --- full posterior (@info$trees_posterior/trees_with_traits_posterior/
     # parameters_posterior), when present. trees_posterior and
     # trees_with_traits_posterior are independent, optional fields, compared
     # the same way via comparePosteriorTreeLists(). ---
     for(prefix in c("trees_posterior", "trees_with_traits_posterior")){
       if(!prefix %in% names(dtreea@info) && !prefix %in% names(dtreeb@info)){
        next
       }
       if(length(dtreea@info[[prefix]]) != length(dtreeb@info[[prefix]])){
         stop(paste(a$clone_id[r], prefix, "not the same length"))
       }
       posterior_comp <- sapply(1:length(dtreea@info[[prefix]]), function(x)
         treesEquivalent(
           dtreea@info[[prefix]][[x]], 
           dtreeb@info[[prefix]][[x]],
           edge_tol))
       if(sum(!posterior_comp) > 0){
         stop(paste(a$clone_id[r], prefix, "not equivalent"))
       }
       treecheck <- treecheck + length(posterior_comp)
     }
  
     ppa <- dtreea@info$parameters_posterior
     ppb <- dtreeb@info$parameters_posterior
     if(is.null(ppa) != is.null(ppb)){
       stop(paste(a$clone_id[r], "parameters_posterior presence not the same"))
     }
     if(!is.null(ppa)){
       if(nrow(ppa) != nrow(ppb)){
         stop(paste(a$clone_id[r], "parameters_posterior row count not the same"))
       }
  
       # Sample can round-trip as "integer" on one side but "double" on the other
       key_a <- paste(sprintf("%.0f", ppa$Sample), ppa$parameter)
       key_b <- paste(sprintf("%.0f", ppb$Sample), ppb$parameter)
       if(!setequal(key_a, key_b)){
         stop(paste(a$clone_id[r], "parameters_posterior Sample/parameter keys not the same"))
       }
       ppb <- ppb[match(key_a, key_b),]
       ok <- isTRUE(all.equal(ppa$value, ppb$value)) #||
         #all(is.na(ppa$value) == is.na(ppb$value) &
         #  (is.na(ppa$value) | ppa$value == ppb$value))
       if(!ok){
         stop(paste(a$clone_id[r], "parameters_posterior values not the same"))
       }else{
         treecheck <- treecheck + nrow(ppa)
       }
     }
    }
    treecheck
  #}
  },error=function(e)e))

  errors <- sapply(checks, function(x) inherits(x, "error"))
  if(sum(errors) == 0){
    treecheck <- sum(unlist(checks))
  }else{
    print(paste(checks[errors]))
    stop("Objects not equivalent")
  }
  if(verbose){
    print(paste("Objects equivalent:",treecheck,"tree checks, 0 failures"))
  }
  return(0)
}
