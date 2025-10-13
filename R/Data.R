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
#' A small example database subset from Turner, J. S. et al. Human germinal centres 
#' engage memory and naive B cells after influenza vaccination. Nature 586, 127–132 (2020).
#'
#' @format   A data.frame with the following Change-O style columns:
#'   \itemize{
#'     \item  \code{sequence_id}:           Sequence identifier
#'     \item \code{sequence}:               B cell sequence
#'     \item \code{productive}:             A logical indicating if the sequence is productive.
#'     \item  \code{v_call}:                V region allele assignments.
#'     \item \code{d_call}:                 D region allele assignments. 
#'     \item  \code{j_call}:                J region allele assignments.
#'     \item  \code{sequence_alignment}:    Sequence alignment.
#'     \item \code{germline_alignment}:     Germline alignment without gaps. 
#'     \item \code{junction}:               Junction
#'     \item \code{juncation_aa}:           Junction aa
#'     \item \code{vj_inframe}:             A logical to see if the vj genes are in frame
#'     \item \code{stop_codon}:             A indicator if there is a stop codon within the alignment
#'     \item \code{locus}:                  Locus identifier. 
#'     \item \code{v_sequence_start}:       Where the V gene starts
#'     \item \code{v_sequence_end}:         Where the V gene ends
#'     \item \code{v_germline_start}:       Where the V germline starts
#'     \item \code{v_germline_end}:         Where the V germline ends
#'     \item \code{np1_length}:             Length of np1
#'     \item \code{d_sequence_start}:       Where the D gene starts
#'     \item \code{d_sequence_end}:         Where the D gene ends
#'     \item \code{d_germline_start}:       Where the D germline starts
#'     \item \code{d_germline_end}:         Where the D germline ends
#'     \item \code{np2_length}:             Length of np2
#'     \item \code{j_sequence_start}:       Where the J gene starts
#'     \item \code{j_sequence_end}:         Where the J gene ends
#'     \item \code{j_germline_start}:       Where the J germline starts
#'     \item \code{j_germline_end}:         Where the J germline ends
#'     \item  \code{junction_length}:       Length of the junction region in nucleotides.
#'     \item \code{v_score}:                V score
#'     \item \code{v_identity}:             Identity score of V
#'     \item \code{v_support}:              V support
#'     \item \code{d_score}:                D score
#'     \item \code{d_identity}:             D identity 
#'     \item \code{d_support}:              D support
#'     \item \code{j_score}:                J score
#'     \item \code{j_support}:              J support
#'     \item \code{j_identity}:             J identity 
#'     \item \code{cell_id}:                Cell identifier 
#'     \item \code{consensus_count}:        Consensus count 
#'     \item \code{indels}:                 Logical if indels are present 
#'     \item \code{sequence_vdj}:           VDJ sequence
#'     \item \code{v_germ_start_vdj}:       Where the V germline starts on the VDJ
#'     \item \code{v_germ_end_vdj}:         Where the V germline ends on the VDJ
#'     \item \code{subject}:                Subject identifier 
#'     \item \code{timepoint}:              Day the sample was taken 
#'     \item \code{cell_type}:              Type of cell 
#'     \item \code{replicate}:              Replicate number 
#'     \item  \code{clone_id}:              Change-O assignment clonal group identifier.
#'     \item \code{seq_type}:               Identifier of data type (10x)
#'     \item \code{vj_gene}:                VJ gene
#'     \item \code{vj_alt_gene}:            Alternative VJ gene
#'     \item \code{v_germline_length}:      Length of the V germline segment
#'     \item \code{d_germline_length}:      Length of the D germline segment 
#'     \item \code{j_germline_lenght}:      Length of the J germline segment 
#'     \item \code{germline_alignment_d_mask}:  Germline alignment with gaps
#' }
"ExampleMixedDb"

#' Example Multiple Partition Trees
#'
#' A small example database subset from Turner, J. S. et al. Human germinal centres 
#' engage memory and naive B cells after influenza vaccination. Nature 586, 127–132 (2020).
#'
#' @format   A data.frame with the following Change-O style columns:
#'   \itemize{
#'     \item  \code{clone_id}:           Clonal cluster
#'     \item  \code{data}:               List of airrClone objects
#'     \item  \code{locus}:              Locus identifier.
#'     \item  \code{seqs}:               Number of sequences
#'     \item \code{igphyml_partitioned_trees}:      IgPhyML partitioned tree
#'     \item \code{raxml_partitioned_trees}:        RAxML partitioned tree
#' }
"ExampleMixedClones"

#' Example AIRR database for TyCHE
#'
#' A small example database from simble.
#'
#' @format   A data.frame with the following AIRR style columns:
#'   \itemize{
#'     \item  \code{sequence_id}:           Sequence identifier
#'     \item  \code{sequence_alignment}:    IMGT-gapped observed sequence.
#'     \item  \code{germline_alignment}:    IMGT-gapped germline sequence.
#'     \item  \code{v_call}:                V region allele assignments.
#'     \item  \code{d_call}:                D region allele assignments.
#'     \item  \code{j_call}:                J region allele assignments.
#'     \item  \code{junction}:              Junction region sequence.
#'     \item  \code{junction_length}:       Length of the junction region in nucleotides.
#'     \item  \code{np1_length}:            Combined length of the N and P regions proximal
#'                                          to the V region.
#'     \item  \code{np2_length}:            Combined length of the N and P regions proximal
#'                                          to the J region.
#'     \item  \code{sample_time}:           Time point of the sample.
#'     \item  \code{location}:              Location of tissue from which the sample was taken.
#'     \item  \code{clone_id}:              Clonal group identifier.
#' }
#' @references
#' \enumerate{
#'   \item  Pre-submission.
#' }
"ExampleAirrTyCHE"

NULL
