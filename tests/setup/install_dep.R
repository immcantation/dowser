#!/usr/bin/env Rscript

# Install dependencies from CRAN or Bitbucket as needed.

####
# Add here Biocondutor dependencies that 
# are not installed in the immcantation/test container
bioconductor_deps <- NULL
####

library(devtools)
library(versions)
# All immcantation R packages
immcantation_packages <- c("alakazam", "shazam", "tigger", "scoper", "dowser", "enchantr")

# Function to install dependencies
installDep <- function(pkg, devel_mode, immcantation=immcantation_packages, 
                       repos="http://lib.stat.cmu.edu/R/CRAN/") {
    
    # Required version 
    pkg_name <- strsplit(pkg," ")[[1]][1]
    pkg_version <- gsub(".*\\([^0-9.]*(.*)\\)$", "\\1", pkg)
    pkg_logic <- gsub("(.*)( *)\\(([><=]*)[^0-9.]*(.*)\\)$", "\\3", pkg)
    if (pkg_version == pkg_name ) { 
        pkg_version <- NULL
        pkg_logic <- ">="
    }
    
    is_immcantation <- pkg_name %in% immcantation
    
    message("Package requested: ", pkg)
    message("devel_mode: ", devel_mode)
    message("is_immcantation: ", is_immcantation)
    
    # Installed packages
    is_installed <- F
    installed_version <- NA
    installed_packages <- installed.packages()
    installed_pkg <-  installed_packages[installed_packages[,"Package"]==pkg_name,,drop=F]
    
    if (nrow(installed_pkg) > 0 ) {
        # Check version of installed package
        tmp_pkg_version <- pkg_version
        if (is.null(pkg_version)) {
            tmp_pkg_version <- "0.0.0"
        }
        installed_version <- numeric_version(installed_pkg[,"Version"])
        is_installed <- eval(parse(text=paste0("numeric_version('",installed_version,"') ",
                                               pkg_logic,
                                               " numeric_version('",tmp_pkg_version,"')")))
    }
    
    if (!is_immcantation | !devel_mode) {
        if (is_installed) {
            message(pkg, " is available.")
        } else {
            # Install from CRAN
            tryCatch({ devtools::install_version(pkg_name, pkg_version, repos="http://lib.stat.cmu.edu/R/CRAN/") },
                     error=function(e) { 
                         # This is needed if there is an Immcantation release package that is not 
                         # available from CRAN
                         cat(e, "\n")
                         message("Installing ",pkg," from Bitbucket...\n ")
                         install_bitbucket(paste0("kleinstein/", pkg_name, "@",pkg_version))
                     })
        }
    } else {
        message(paste0(pkg,": installing most recent version from Bitbucket @master.")) 
        install_bitbucket(paste0("kleinstein/", pkg_name, "@master"), upgrade = "never")
    }
}

# Parse this package version in DESCRIPTION
this_pkg_version <- read.dcf("DESCRIPTION", fields="Version")
this_pkg_version <- gsub(".*\\([^0-9.]*(.*)\\)$", "\\1", this_pkg_version)

# If the package is using the devel version number (ends in .999)
# always install devel versions of immcantation packages
devel_mode <- grepl("\\.999$", this_pkg_version)

# Parse Imports field in DESCRIPTION
d <- read.dcf("DESCRIPTION", fields="Imports")
d <- sub("^\\n", "", d)
imports <- strsplit(d, ",\n")[[1]]

# Install immcantation packages first
immcantation <- unlist(sapply(immcantation_packages, grep, imports))
imports <- imports[unique(c(immcantation, 1:length(imports)))]
# Skip bioconductor packages
imports <- imports[!imports %in% bioconductor_deps]


# Install
for (i in 1:length(imports)) {
    this_import <- imports[i]
    installDep(this_import, devel_mode)
}

if (!is.null(bioconductor_deps)) {
    BiocManager::install(bioconductor_deps)
}