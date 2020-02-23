#### Lineage classes ####

#' S4 class defining a clone
#' 
#' \code{ChangeoClone} defines a common data structure for perform lineage recontruction
#' from Change-O data.
#' 
#' @slot     data        data.frame containing sequences and annotations. Contains the
#'                       columns \code{SEQUENCE_ID} and \code{SEQUENCE}, as well as any additional 
#'                       sequence-specific annotation columns.
#' @slot     clone       string defining the clone identifier.
#' @slot     germline    string containing the heavy chain germline sequence for the clone.
#' @slot     lgermline   string containing the light chain germline sequence for the clone.
#' @slot     hlgermline  string containing the combined germline sequence for the clone.
#' @slot     v_gene      string defining the V segment gene call.
#' @slot     j_gene      string defining the J segment gene call.
#' @slot     junc_len    numeric junction length (nucleotide count).
#' @slot     locus       the loci represented in this clone object
#' @slot     region      index showing which regions represent which site
#' @slot     phylo_seq   sequence column used for phylogenetic tree building
#' @seealso  See \link{makeChangeoClone} and \link{buildPhylipLineage} for use.
#'           
#' @name         ChangeoClone-class
#' @rdname       ChangeoClone-class
#' @aliases      ChangeoClone
#' @exportClass  ChangeoClone
setClass("ChangeoClone", 
         slots=c(data="data.frame",
         clone="character",
         germline="character", 
         lgermline="character",
         hlgermline="character",
         v_gene="character", 
         j_gene="character", 
         junc_len="numeric",
         locus="character",
    	 region="character",
       phylo_seq="character"))

