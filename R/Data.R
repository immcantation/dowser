# Documentation and definitions for data and constants

#### Sysdata ####

# 1x20 vector of default amino acid hydropathy scores
# HYDROPATHY_KYTJ82

# 1x20 vector of default amino acid bulkiness scores
# BULKINESS_ZIMJ68

# 1x20 vector of default amino acid polarity scores
# POLARITY_GRAR74

# 1x7 vector of default amino acid pK values
# PK_EMBOSS

#### Data ####

#' Example AIRR database
#'
#' A small example database subset from Laserson and Vigneault et al, 2014.
#'
#' @format   A data.frame with the following AIRR style columns:
#'   \itemize{
#'     \item  \code{sequence_id}:           Sequence identifier
#'     \item  \code{sequence_alignment}:    IMGT-gapped observed sequence.
#'     \item  \code{germline_alignment_d_mask}:  IMGT-gapped germline sequence with N, P and 
#'                                          D regions masked.
#'     \item  \code{v_call}:                V region allele assignments.
#'     \item  \code{v_call_genotyped}:      TIgGER corrected V region allele assignment.
#'     \item  \code{d_call}:                D region allele assignments.
#'     \item  \code{j_call}:                J region allele assignments.
#'     \item  \code{junction}:              Junction region sequence.
#'     \item  \code{junction_length}:       Length of the junction region in nucleotides.
#'     \item  \code{np1_length}:            Combined length of the N and P regions proximal
#'                                          to the V region.
#'     \item  \code{np2_length}:            Combined length of the N and P regions proximal
#'                                          to the J region.
#'     \item  \code{sample}:                Sample identifier. Time in relation to vaccination.
#'     \item  \code{isotype}:               Isotype assignment.
#'     \item  \code{duplicate_count}:       Copy count (number of duplicates) of the sequence.
#'     \item  \code{clone_id}:              Change-O assignment clonal group identifier.
#' }
#' 
#' @seealso \link{ExampleDbChangeo} \link{ExampleTrees}
#' 
#' @references
#' \enumerate{
#'   \item  Laserson U and Vigneault F, et al. High-resolution antibody dynamics of 
#'            vaccine-induced immune responses. 
#'            Proc Natl Acad Sci USA. 2014 111:4928-33.
#' }
"ExampleDb"

#' Example Change-O database
#'
#' A small example database subset from Laserson and Vigneault et al, 2014.
#'
#' @format   A data.frame with the following Change-O style columns:
#'   \itemize{
#'     \item  \code{SEQUENCE_ID}:           Sequence identifier
#'     \item  \code{SEQUENCE_IMGT}:         IMGT-gapped observed sequence.
#'     \item  \code{GERMLINE_IMGT_D_MASK}:  IMGT-gapped germline sequence with N, P and 
#'                                          D regions masked.
#'     \item  \code{V_CALL}:                V region allele assignments.
#'     \item  \code{V_CALL_GENOTYPED}:      TIgGER corrected V region allele assignment.
#'     \item  \code{D_CALL}:                D region allele assignments.
#'     \item  \code{J_CALL}:                J region allele assignments.
#'     \item  \code{JUNCTION}:              Junction region sequence.
#'     \item  \code{JUNCTION_LENGTH}:       Length of the junction region in nucleotides.
#'     \item  \code{NP1_LENGTH}:            Combined length of the N and P regions proximal
#'                                          to the V region.
#'     \item  \code{NP2_LENGTH}:            Combined length of the N and P regions proximal
#'                                          to the J region.
#'     \item  \code{SAMPLE}:                Sample identifier. Time in relation to vaccination.
#'     \item  \code{ISOTYPE}:               Isotype assignment.
#'     \item  \code{DUPCOUNT}:              Copy count (number of duplicates) of the sequence.
#'     \item  \code{CLONE}:                 Change-O assignment clonal group identifier.
#' }
#' 
#' @seealso \link{ExampleDb} \link{ExampleTrees}
#' 
#' @references
#' \enumerate{
#'   \item  Laserson U and Vigneault F, et al. High-resolution antibody dynamics of 
#'            vaccine-induced immune responses. 
#'            Proc Natl Acad Sci USA. 2014 111:4928-33.
#' }
"ExampleDbChangeo"


#' Example Ig lineage trees
#'
#' A set of Ig lineage trees generated from the \code{ExampleDb} file, subset to
#' only those trees with at least four nodes.
#'
#' @format   A list of igraph objects output by alakazam::buildPhylipLineage.
#'           Each node of each tree has the following annotations (vertex attributes):
#'   \itemize{
#'     \item  \code{sample}:    Sample identifier(s). Time in relation to vaccination.
#'     \item  \code{isotype}:   Isotype assignment(s). 
#'     \item  \code{duplication_count}:  Copy count (number of duplicates) of the sequence.
#'   }
#'   
#' @seealso \link{ExampleTrees}
"ExampleTrees"


NULL
