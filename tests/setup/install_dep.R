#!/usr/bin/env Rscript

# Install dependencies from CRAN or Bitbucket as needed.

# Imports
library(devtools)
library(versions)

# Function to install packages
installDep <- function(this_pack_v, dep_pack_name, dep_pack_v) {
    required_version <- gsub(".*\\([^0-9.]*(.*)\\)$", "\\1", dep_pack_v)
    devel <- length(grep("\\.999$", required_version)) > 0
    
    get_cran_versions <- function (pkg) {
        
        current_url <- sprintf("%s/src/contrib",versions:::latest.MRAN())
        lines <- versions:::url.lines(current_url)
        lines <- lines[grep("^<a href=\"*", lines)]
        tarballs <- gsub(".*href=\"([^\"]+)\".*", "\\1", lines)
        dates <- gsub(".*  ([0-9]+-[a-zA-Z]+-[0-9]+) .*", "\\1", lines)
        dates <- as.Date(dates, format = "%d-%b-%Y")
        idx <- grep(sprintf("^%s_.*.tar.gz$", pkg), tarballs)
        # Keep all versions
        if (length(idx) < 1) {
            warning(sprintf("The current version and publication date of %s could not\n                     be detected", 
                            pkg))
            versions <- dates <- NA
        } else if (length(idx) > 0) {
            versions <- tarballs[idx]
            versions <- gsub(sprintf("^%s_", pkg), "", versions)
            versions <- numeric_version(gsub(".tar.gz$", "", versions))
            dates <- dates[idx]
        } 
        
        ret <- list()
        ret[[pkg]] <- data.frame(
                version = versions, 
                date = as.character(dates), 
                stringsAsFactors = FALSE)
        ret
    }
    
    cran_versions <- get_cran_versions(dep_pack_name)
    cran_versions <- cran_versions[[dep_pack_name]]$version
    
    this_pack_devel <- length(grep("\\.999$", this_pack_v)) > 0
    in_cran <- numeric_version(required_version) %in% cran_versions
    
    if (!this_pack_devel & !devel & in_cran) {
        tryCatch({ devtools::install_version(dep_pack_name, required_version, repos="http://lib.stat.cmu.edu/R/CRAN/") },
                 error=function(e) { 
                     cat(e, "\n")
                     message("Installing from Bitbucket...\n ")
                     install_bitbucket(paste0("kleinstein/", dep_pack_name, "@master"))
                 })
    } else {
        if (!in_cran & !devel) { 
            warning(paste0(required_version," not found in CRAN.")) 
        }
        message(paste0(dep_pack_name, " ", required_version,": installing most recent version from Bitbucket.")) 
        install_bitbucket(paste0("kleinstein/", dep_pack_name, "@master"), upgrade = "never")
    }
}

# Parse Imports field in DESCRIPTION
pkg_version <- read.dcf("DESCRIPTION", fields="Version")
d <- read.dcf("DESCRIPTION", fields="Imports")
d <- sub("^\\n", "", d)
imports <- strsplit(d, ",\n")[[1]]

# Install
idx <- sapply(c("shazam", "alakazam"), grep, imports)
for (i in 1:length(idx)) {
    this_package_name <-  names(idx)[[i]]
    this_package_version <-  imports[idx[[i]]]
    installDep(pkg_version, this_package_name, this_package_version)
}
