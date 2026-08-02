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
#' @param   dowser_fields include dowser-specific information? (recommended)
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
  edge_tol=1e-8, cell="cell_id", heavy="IGH", dowser_fields=TRUE){
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
      clone$program_origin <- "dowser"
      #clone$info$locus <- object$data[[row]]@locus
      clone$info$region <- object$data[[row]]@region
      clone$info$numbers <- object$data[[row]]@numbers
      clone$info$phylo_seq <- object$data[[row]]@phylo_seq
      clone$v_call <- object$data[[row]]@v_gene
      clone$j_call <- object$data[[row]]@j_gene
      clone$junction_length <- object$data[[row]]@junc_len
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

    phy <- object$trees[[row]]
    phy$node.label <- paste0("Node",(length(phy$tip.label)+1):length(phy$nodes),"-",clone$clone_id)

    phy$tip.label[phy$tip.label == "Germline"] = paste0("Germline-",clone$clone_id)

    # add info from other columns
    forbidden <- c("region", "numbers", "phylo_seq", "germline",
      "lgermline", "hlgermline","program_origin")
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
    #clone$info$data <- tipdata

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
    germline_node <- which(phy$tip.label == paste0("Germline-",clone$clone_id))

    # metadata columns for tips
    nodes <- list()
    for(i in 1:length(phy$nodes)){
      node <- list()
      if(dowser_fields){
        node$info <- list()
      }
      node$node_class <- node_class
      node$repertoire_id <- repertoire_id
      receptors <- list()
      if(i <= length(phy$tip.label) && i != germline_node){
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
        if(i == germline_node){
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

      seqa <- getNodeSeq(object, node=i, tree=phy, gaps=TRUE)
      seq <- getNodeSeq(object, node=i, tree=phy, gaps=FALSE)
      ga <- getNodeSeq(object, node=germline_node, tree=phy, gaps=TRUE)
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

      if(node$number == germline_node){
        clone$inferred_ancestor <- node$node_id
      }
    }
    clone$nodes <- nodes

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
    nobject <- readTreesJSON(file)
    # dowserObjectEquivalent understands pml fits directly (reducing them via
    # pmlParamsToList itself), so no pre-normalization is needed here.
    validate <- dowserObjectEquivalent(object, nobject, verbose, edge_tol, dowser_fields)
  }
}


#' \code{readTreesJSON}
#' Experimental. Read trees from JSON/AIRR format from Dowser
#' @param    file   .json file
#' @param    heavy    name of heavy chain locus
#' @param    light    names of light chain loci
#' @details
#' Reads files written by \code{\link{writeTreesJSON}}, including trees built
#' with any \code{getTrees} \code{build} option.
#' @export
readTreesJSON = function(file, heavy="IGH", light=c("IGK","IGL")){
  rclones <- jsonlite::read_json(file)
  output <- dplyr::tibble()
  outtrees <- list()
  outdata <- list()
  rearrangements <- rclones$Rearrangement

  if("program_origin" %in% names(rclones$Clone[[1]]) &&
    rclones$Clone[[1]]$program_origin == "dowser"){
    dowser <- TRUE
  }else{
    dowser <- FALSE
    warning("Some features limited by using non-Dowser JSON")
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
    ancestor_node_id <- clone$inferred_ancestor

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
        }else if(heavy %in% loci && length(intersect(light,unique_loci)) > 0){
          phylo_seq <- "hlsequence"
        }else if(unique_loci[1] == heavy && length(unique_loci) == 1){
          phylo_seq <- "sequence"
        }else if(length(intersect(light,unique_loci)) > 0){
          phylo_seq <- "lsequence"
        }else{
          stop(paste("Could not assign phylo_seq based on locus column:",
            paste(unique_loci, collapse=",")))
        }
      }
      if(nchar(seq) != length(loci)){
        stop("seq and loci vector not the same length!")
      }
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
          if(dowser){ #only mess with the germline name if coming from Dowser
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
    for(i in 1:length(tips)){
      tip <- tips[[i]]
      if(dowser){
        row <- tip$info$tipdata
      }else{
        row <- list()
        row$sequence_id <- tip$node_id
      }
      row[[phylo_seq]]<- tip$sequence

      if(tip$node_id != ancestor_node_id){
        data[[i]] <- row
      }else{
        germline <- row[[phylo_seq]]
      }
    }
    bdata <- dplyr::bind_rows(data)

    if(dowser){
      bdata <- bdata[order(bdata$tip_order),]
      bdata <- dplyr::select(bdata, -"tip_order")
      rphy$tip.label[rphy$tip.label == ancestor_node_id] <- "Germline"
    }

    lgermline <- ""
    hlgermline <- ""
    if(!dowser){
      region <- rep("N",times=nchar(germline))
    }
    if(phylo_seq == "hlsequence"){
      hlgermline <- germline
      lgermline <- ""
      germline <- ""
      if(!dowser){
        region <- rep("N",times=nchar(hlgermline))
      }
    }else if(phylo_seq == "lsequence"){
      lgermline <- germline
      germline <- ""
      hlgermline <- ""
      if(!dowser){
        region <- rep("N",times=nchar(lgermline))
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
    outtrees[[ci]] <- rphy

    # store info in other columns
    temp <- dplyr::tibble(clone_id=clone$clone_id[[1]], seqs=nrow(outclone@data),
      locus=paste0(sort(unique(loci)),collapse=","))
    for(n in names(clone$info)){
      if(!n %in% c("region", "numbers", "phylo_seq", "germline", "lgermline",
        "hlgermline", "trees", "data", "clone_id", "seqs", "locus", "program_origin", "treeinfo")){
        ni <- clone$info[[n]]
        if(length(ni) > 1){
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
  standard <- c("clone_id", "data", "seqs", "locus", "trees")
  nonstandard <- setdiff(onames, standard)

  output <- dplyr::select(output, dplyr::any_of(standard), dplyr::any_of(nonstandard))
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
      if(!all.equal(va, vb, tolerance=tol)){
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
    if(!all.equal(ca, cb, tolerance=edge_tol, check.attributes=FALSE)){
      return(FALSE)
    }
  }

  return(TRUE)
}

#'\code{dowserObjectEquivalent}
#' Experimental. Check if two Dowser objects are equivalent
#' @param    obj1  First Dowser object
#' @param    obj2  Second Dowser object
#' @param    verbose  print out more info
#' @param    edge_tol tolerance for branch length checks (if check=TRUE)
#' @param   dowser_fields check dowser-specific fields and gapped sequences?
#' @details
#' In addition to the existing tree topology, edge length, sequence, and
#' data slot checks, this also verifies \code{tree$parameters} when present
#' -- including \code{build="pml"} trees.
#' @export
dowserObjectEquivalent = function(obj1, obj2, verbose=TRUE, edge_tol=1e-8,
  dowser_fields=TRUE){
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

  treecheck <- 0
  datacheck <- 0
  nodewarn <- 0
  for(row in 1:nrow(a)){
    # check trees
    treea <- a$trees[[row]]
    treeb <- b$trees[[row]]

    tipsa <- treea$tip.label
    tipsb <- treeb$tip.label
    if(!setequal(tipsa, tipsb)){
      diffs <- c(setdiff(tipsb, tipsa), setdiff(tipsa, tipsb))
      fdiffs <- diffs[!diffs %in% c("Germline",paste0("Germline-",a$clone_id))]
      if(!dowser_fields && length(fdiffs) == 0){
        treea$tip.label[treea$tip.label %in% diffs] = "Germline"
        treeb$tip.label[treeb$tip.label %in% diffs] = "Germline"
      }else{
        stop(paste(a$clone_id[row],"tips not the same"))
      }
    }else{
      treecheck <- treecheck + 1
    }
    all_subtrees_a <- lapply(1:length(treea$nodes), function(x)getSubTaxa(x, treea))
    all_subtrees_b <- lapply(1:length(treeb$nodes), function(x)getSubTaxa(x, treeb))
    if(length(all_subtrees_a) != length(all_subtrees_b)){
      stop(paste(a$clone_id[row],"trees don't have the same subtree numbers"))
      treecheck <- treecheck + 1
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
        stop(paste(a$clone_id[row],"Subtrees don't match"))
      }else{
        treecheck <- treecheck + 1
      }
      seqa <- getNodeSeq(a, clone=a$clone_id[row], node=sta, gaps=dowser_fields)
      seqb <- getNodeSeq(b, clone=b$clone_id[row], node=match, gaps=dowser_fields)
      if(!is.null(seqa)){
        if(is.null(seqb)){
          stop(paste(a$clone_id[row], "One seq is NULL", sta, match))
        }
        if(length(seqa) != length(seqb)){
          stop(paste(a$clone_id[row], "different seq numbers", sta, match))
        }
        for(seqindex in 1:length(seqa)){
          if(!seqa[seqindex] == seqb[seqindex]){
            stop(paste(a$clone_id[row], "Seqs don't match", sta, match))
          }else{
            treecheck <- treecheck + 1
          }
        }
      }
      ea <- round(treea$edge.length[treea$edge[,2] == sta], digits=11)
      eb <- round(treeb$edge.length[treeb$edge[,2] == match], digits=11)
      if(length(ea) > 0 || length(eb) > 0){
        if(abs(ea - eb) > edge_tol){
          stop(paste(a$clone_id[row], "Edges not within edge_tol",ea, eb, sta, match))
        }else{
          treecheck <- treecheck + 1
        }
      }
      #to do: add check for node names
      if(!is.null(treea$state) || !is.null(treeb$state)){
        if(treea$state[sta] != treeb$state[match]){
          stop(paste(a$clone_id[row], "node states not equal",treea$state[sta],
            treea$state[match], sta, match))
        }
      }
    }
    # check remaining tree info
    if(dowser_fields && !setequal(names(treea), names(treeb))){
      namediff <- setdiff(names(treea), names(treeb))
      if(length(namediff) > 1 || namediff != "node.label"){
        print(namediff)
        stop(paste(a$clone_id[row], "Tree names not equal"))
      }
      if("node.label" %in% c(names(treea),names(treeb))){
        if(nodewarn == 0){
          warning("node.label not currently checked, or preserved by writeTreesJSON")
        }
        nodewarn <- nodewarn + 1
      }
    }else{
      treecheck <- treecheck + 1
    }
    for(n in names(treea)){
      if(n %in% c("edge","tip.label","Nnode","edge.length","nodes","node.label")){
        next
      }
      if(n == "parameters"){
        if(is.null(treea[[n]]) && is.null(treeb[[n]])){
          treecheck <- treecheck + 1
        }else if(is.null(treea[[n]]) || is.null(treeb[[n]])){
          stop(paste(a$clone_id[row], n, "not the same"))
        }else if(!pmlParamsEqual(treea[[n]], treeb[[n]], edge_tol=edge_tol)){
          stop(paste(a$clone_id[row], n, "not the same"))
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
        stop(paste(a$clone_id[row], n, "not the same"))
      }else{
        treecheck <- treecheck + 1
      }
    }

    # check data
    da <- a$data[[row]]
    db <- b$data[[row]]
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
          if(sum(dataa[[name]] != datab[[name]]) != 0){
            stop(paste(a$clone_id[row], "data", name, "not the same"))
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
        stop(paste(a$clone_id[row], n, "slot is not equal"))
      }else{
        datacheck <- datacheck + length(slot(da,n))
      }
    }
  }
  if(verbose){
    print(paste("Objects equivalent:",treecheck,"tree checks,",datacheck,"data checks, 0 failures"))
  }
  return(0)
}
