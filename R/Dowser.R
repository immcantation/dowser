# Dowser package documentation and import directives

#' The dowser package
#' 
#' \code{dowser} is a phylogenetic analysis package as part of the Immcantation suite of tools.
#' For additional details regarding the use of the \code{dowser} package see the 
#' vignettes:\cr
#' \code{browseVignettes("dowser")}
#' 
#' @name     dowser
#' @docType  package
#' @references
#' \enumerate{
#'   \item  Hoehn, KB, Pybus, OG, Kleinstein SH (submitted) Phylogenetic analysis of 
#' 		migration, differentiation, and class switching	in B cells.

#' }
#' 
#' @import      ggplot2
#' @import      ggtree
#' @import      graphics
#' @import      methods
#' @import      utils
#' @importFrom  tidyselect  all_of
#' @importFrom  gridExtra   grid.arrange
#' @importFrom  stringr     str_count    
#' @importFrom  alakazam    makeChangeoClone DNA_IUPAC
#' @importFrom  phylotate   read_annotated
#' @importFrom	ape 		read.tree di2multi reorder.phylo root ladderize subtrees
#' @importFrom  dplyr       do n desc %>% pull
#'                          bind_cols bind_rows combine arrange 
#'                          group_by ungroup
#'                          filter slice select 
#'                          mutate mutate_at 
#' 							one_of if_else
#'							right_join rowwise
#'                          summarize summarize_at
#'                          transmute rename

NULL
