library(markr)
library(dowser)

# Directories
pkg_path <- "."
doc_path <- "./docs"

# Build
build_mkdocs(pkg_path, doc_path=doc_path, yaml=F)
