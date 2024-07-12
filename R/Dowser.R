#' @keywords internal
#' @aliases dowser-package
"_PACKAGE"


# Dowser package documentation and import directives
#'
#' The dowser package
#' 
#' \code{dowser} is a phylogenetic analysis package as part of the Immcantation suite of tools.
#' For additional details regarding the use of the \code{dowser} package see the 
#' vignettes:\cr
#' \code{browseVignettes("dowser")}
#' 
#' @name     dowser
#' @references
#' \enumerate{
#'   \item  Hoehn KB, Pybus OG, Kleinstein SH (2022) Phylogenetic analysis of 
#'         migration, differentiation, and class switching in B cells.
#'         PLoS Computational Biology. https://doi.org/10.1371/journal.pcbi.1009885
#' }
#' 
#' @import      ggplot2
#' @import      ggtree
#' @import      graphics
#' @import      markdown
#' @import      methods
#' @import      utils
#' @importFrom  shazam      setRegionBoundaries
#' @importFrom  phangorn    phyDat acctran pratchet ancestral.pars dist.ml
#'                          NJ pml optim.pml pml.control ancestral.pml
#' @importFrom  tidyselect  all_of
#' @importFrom  gridExtra   grid.arrange
#' @importFrom  stringr     str_count    
#' @importFrom  alakazam    checkColumns maskSeqGaps padSeqEnds 
#'                          collapseDuplicates getGene translateDNA seqDist
#'                          getAllele buildPhylipLineage graphToPhylo DNA_IUPAC
#'                          makeTempDir
#' @importFrom  phylotate   read_annotated
#' @importFrom  ape         read.tree di2multi reorder.phylo root ladderize subtrees
#'                          as.AAbin as.DNAbin getMRCA dist.nodes multi2di extract.clade
#'                          keep.tip bind.tree collapse.singles bind.tree unroot read.dna
#'                          is.rooted
#' @importFrom  dplyr       do n desc %>% pull
#'                          bind_cols bind_rows combine arrange 
#'                          group_by ungroup
#'                          filter slice select 
#'                          mutate mutate_at 
#'                          one_of
#'                          right_join rowwise
#'                          summarize summarize_at
#'                          transmute rename
#' @importFrom Biostrings   pairwiseAlignment
#' @importFrom airr		    read_rearrangement

NULL

"_PACKAGE"