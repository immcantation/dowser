# Download and Installation

Download
-------------------------------------------------------------------------------

The latest stable release of `dowser` can be downloaded from 
<a href="http://cran.rstudio.com/web/packages/dowser" target="_blank">CRAN</a>
or <a href="https://github.com/immcantation/dowser/tags" target="_blank">GitHub</a>.

Installing Released Versions
-------------------------------------------------------------------------------

The simplest method is to use Bioconductor's `install` function, which will install Bioconductor dependencies:

```R

install.packages("BiocManager")

BiocManager::install("dowser")

```

You can also install `dowser` via CRAN, but Bioconductor dependencies may not install:

```R
install.packages("dowser")
```

Downloaded source builds from GitHub may be installed in the usual way:

```R
install.packages("dowser_x.y.z.tar.gz", repos = NULL, type = "source")
```

If you have any trouble installing the package, it may be due to the Bioconductor 
dependencies. You can run the following command to see what other packages may be needed:

```R
available.packages()["dowser", "Imports"]
```

Try checking error messages and installing any failed dependencies individually. 
You may need to upgrade to the latest version of R. If you have questions, you 
can email the [Immcantation Group](mailto:immcantation@googlegroups.com).

Building Development Versions
-------------------------------------------------------------------------------

To build from the [source code](http://github.com/immcantation/dowser),
first install the build dependencies:

```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "Rcpp"))
```

To install the latest development code via devtools:

```R
library(devtools)
install_github("immcantation/dowser@master")
```

Note, using `install_github` will not build the documentation. To generate the 
documentation, clone the repository and build as normal using `devtools`, 
`roxygen` and `knitr`:

```R
library(devtools)
install_deps()
document()
build()
install()
```
