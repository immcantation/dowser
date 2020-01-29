# Dowser package documentation and import directives

#' The dowser package
#' 
#' \code{dowser} in a member of the Immcantation suite of tools and serves five main 
#' purposes:
#' \itemize{
#'   \item  Providing core functionality for other R packages in the Change-O suite. This
#'          includes common tasks such as file I/O, basic DNA sequence manipulation, and
#'          interacting with V(D)J segment and gene annotations.
#'   \item  Providing an R interface for interacting with the output of the pRESTO 
#'          tool suite.
#'   \item  Performing lineage reconstruction on clonal populations of immunoglobulin 
#'          (Ig) sequences. 
#'   \item  Performing clonal abundance and diversity analysis on lymphocyte repertoires.
#'   \item  Performing physicochemical property analyses of lymphocyte receptor sequences.
#' }
#' For additional details regarding the use of the \code{dowser} package see the 
#' vignettes:\cr
#' \code{browseVignettes("dowser")}
#' 
#' @name     dowser
#' @docType  package
#' @references
#' \enumerate{
#'   \item  Hoehn, KB, Pybus, OG, Kleinstein SH (submitted) Phylogenetic analysis of 
#' 		migration, differentiation, and class switching	in B cells using maximum parsimony.

#' }
#' 
#' @import      ggplot2
#' @import      ggtree
#' @import      graphics
#' @import      methods
#' @import      utils
#' @importFrom  alakazam    makeChangeoClone
#' @importFrom  phylotate   read_annotated
#' @importFrom	ape 		read.tree di2multi reorder.phylo root ladderize
#' @importFrom  dplyr       do n desc %>%
#'                          bind_cols bind_rows combine arrange 
#'                          group_by ungroup
#'                          filter slice select 
#'                          mutate mutate_at 
#' 							one_of if_else
#'							right_join rowwise
#'                          summarize summarize_at
#'                          transmute rename

NULL
