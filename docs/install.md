Download

Building Development Versions
-------------------------------------------------------------------------------

Currently dowser can only be built from the [source code](http://bitbucket.org/kleinstein/dowser),
first install the build dependencies:

```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "Rcpp"))
```

To install the latest development code via devtools, along with a development branch of alakazam:

```R
library(devtools)
install_bitbucket("kleinstein/alakazam@dowser")
install_bitbucket("kleinstein/dowser@master")
```

Note, using `install_bitbucket` will not build the documentation. To generate the 
documentation, clone the repository and build as normal using devtools, 
roxygen and knitr:

```R
library(devtools)
install_deps()
document()
build()
install()
```
