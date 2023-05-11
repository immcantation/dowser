##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  dowser_ref <- paste('Hoehn K, Pybus O, Kleinstein S (2022). Phylogenetic',
                      'analysis of migration, differentiation, and class switching in B cells. PLoS',
                      'Computational Biology. https://doi.org/10.1371/journal.pcbi.1009885.')
  ref <- paste0("If you are using dowser in published research please cite ", dowser_ref)
  if(!is.null(ref)) packageStartupMessage(ref)
}