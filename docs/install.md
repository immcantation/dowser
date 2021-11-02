Download
-------------------------------------------------------------------------------

Dowser can be built from the [source code](http://bitbucket.org/kleinstein/dowser), or used in the 
Immcantation Docker container development build. The Docker container currently the only way to use 
the igphyml dependent features of Dowser on Windows.

Using Docker image
-------------------------------------------------------------------------------

To use Dowser in the Docker image, simply pull and run the immcantation/suite:devel 
[Docker image](https://immcantation.readthedocs.io/en/stable/docker/intro.html). Instructions in 
hyperlink. Note this might not always be the latest version of Dowser.


Building Development Versions
-------------------------------------------------------------------------------

Alternatively, to install the latest version source code, first install the build dependencies:

```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "Rcpp"))
```

To install the latest development code via devtools, along with the development version of Alakazam:

```R
library(devtools)
install_bitbucket("kleinstein/alakazam@master")
install_bitbucket("kleinstein/dowser@master")
```
