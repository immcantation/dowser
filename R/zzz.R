##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  dowser_ref <- paste("Hoehn K, Pybus O, Kleinstein S (2022). “Phylogenetic analysis of migration, differentiation, and class switching in B cells.” PLoS Computational Biology. doi:10.1371/journal.pcbi.1009885 https://doi.org/10.1371/journal.pcbi.1009885, https://doi.org/10.1371/journal.pcbi.1009885.")
  corr_ref <- paste("Hoehn K, Turner J, Miller F, Jiang R, Ellebedy A, Pybus O, Kleinstein S (2021). “Human B cell lineages associated with germinal centers following influenza vaccination are measurably evolving.” eLife. doi:10.7554/eLife.70873 https://doi.org/10.7554/eLife.70873, https://elifesciences.org/articles/70873.")
  igphyml_ref <- paste("Hoehn K, Van der Heiden J, Zhou J, Lunter G, Pybus O, Kleinstein S (2019). “Repertoire-wide phylogenetic models of B cell molecular evolution reveal evolutionary signatures of aging and vaccination.” PNAS. doi:10.1073/pnas.1906020116 https://doi.org/10.1073/pnas.1906020116, https://www.pnas.org/doi/10.1073/pnas.1906020116.")
  ref <- paste0("If you are using dowser in published research please cite ", dowser_ref, "\n\n",
                "If you use the correlationTest function for measurable evolution, please also cite ", corr_ref,
                "\n\n", "If you use IgPhyML for building trees please also cite ", igphyml_ref)
  if(!is.null(ref)) packageStartupMessage(ref)
}
