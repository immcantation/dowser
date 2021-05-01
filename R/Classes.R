# Classes, Generics and Methods

#' S4 class defining a clone in Dowser
#' 
#' \code{airrClone} defines a common data structure for perform lineage recontruction
#' from AIRR data, based heavily on alakazam::ChangeoClone. 
#' 
#' @slot     data        data.frame containing sequences and annotations. Contains the
#'                       columns \code{SEQUENCE_ID} and \code{SEQUENCE}, as well as any additional 
#'                       sequence-specific annotation columns
#' @slot     clone       string defining the clone identifier
#' @slot     germline    string containing the heavy chain germline sequence for the clone
#' @slot     lgermline   string containing the light chain germline sequence for the clone
#' @slot     hlgermline  string containing the combined germline sequence for the clone
#' @slot     v_gene      string defining the V segment gene call
#' @slot     j_gene      string defining the J segment gene call
#' @slot     junc_len    numeric junction length (nucleotide count)
#' @slot     locus       the loci represented in this clone object
#' @slot     chain       index showing which chain represent which site
#' @slot     region      index showing FWR/CDR region for each site
#' @slot     phylo_seq   sequence column used for phylogenetic tree building
#' @slot     numbers     index (usually IMGT) number of each site in \code{phylo_seq}
#' @seealso  See \link{formatClones} for use.
#'           
#' @name         airrClone-class
#' @rdname       airrClone-class
#' @aliases      airrClone
#' @exportClass  airrClone
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
     chain="character",
     region="character",
     phylo_seq="character",
     numbers="numeric"))
