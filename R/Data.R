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
#' @seealso \link{ExampleDbChangeo} \link{ExampleClones}
#' 
#' @references
#' \enumerate{
#'   \item  Laserson U and Vigneault F, et al. High-resolution antibody dynamics of 
#'            vaccine-induced immune responses. 
#'            Proc Natl Acad Sci USA. 2014 111:4928-33.
#' }
"ExampleAirr"

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
#' @seealso \link{ExampleAirr} \link{ExampleClones}
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
#' A tibble of Ig lineage trees generated from the \code{ExampleAirr} file
#'
#' @format   A tibble of airrClone and phylo objects output by getTrees.
#'   \itemize{
#'     \item  \code{clone_id}:   Clonal cluster
#'     \item  \code{data}:       List of airrClone objects
#'     \item  \code{seqs}:       Number of sequences
#'     \item  \code{trees}:      List of phylo objects
#'   }
#'   
#' @seealso \link{ExampleClones}
"ExampleClones"

#' Example Ig lineage trees with biopsy reconstructions.
#'
#' Same as ExampleClones but with biopsies predicted at internal nodes
#'
#' @format   A tibble of airrClone and phylo objects output by getTrees.
#'   \itemize{
#'     \item  \code{clone_id}:   Clonal cluster
#'     \item  \code{data}:       List of airrClone objects
#'     \item  \code{seqs}:       Number of sequences
#'     \item  \code{trees}:      List of phylo objects
#'   }
#'   
#' @seealso \link{BiopsyTrees}
"BiopsyTrees"

#' Example Ig lineage trees with isotype reconstructions.
#'
#' Same as ExampleClones but with isotypes predicted at internal nodes
#'
#' @format   A tibble of airrClone and phylo objects output by getTrees.
#'   \itemize{
#'     \item  \code{clone_id}:   Clonal cluster
#'     \item  \code{data}:       List of airrClone objects
#'     \item  \code{seqs}:       Number of sequences
#'     \item  \code{trees}:      List of phylo objects
#'   }
#'   
#' @seealso \link{IsotypeTrees}
"IsotypeTrees"

#' Example Ig lineage trees sampled over time.
#'
#' Same as ExampleClones but with timepoint as a trait value
#'
#' @format   A tibble of airrClone and phylo objects output by getTrees.
#'   \itemize{
#'     \item  \code{clone_id}:   Clonal cluster
#'     \item  \code{data}:       List of airrClone objects
#'     \item  \code{seqs}:       Number of sequences
#'     \item  \code{trees}:      List of phylo objects
#'   }
#'   
#' @seealso \link{TimeTrees}
"TimeTrees"


#' Example Change-O database
#'
#' A small example database originally created for the Immcantation package 'scoper.
#'
#' @format   A data.frame with the following Change-O style columns:
#'   \itemize{
#'     \item  \code{sequence_id}:           Sequence identifier
#'     \item  \code{subject_id}:            Subject identifier.
#'     \item  \code{v_call}:                V region allele assignments.
#'     \item  \code{j_call}:                J region allele assignments.
#'     \item  \code{sequence_alignment}:    Sequence alignment.
#'     \item \code{locus}:                  Locus identifier. 
#'     \item \code{cell_id}:                Cell identifier 
#'     \item  \code{junction_length}:       Length of the junction region in nucleotides.
#'     \item  \code{clone_id}:              Change-O assignment clonal group identifier.
#'     \item \code{expected_clone_subgroup}: The expected clone subgroup. Used for unit testing. 
#' }
"ExampleMixedDb"


NULL
