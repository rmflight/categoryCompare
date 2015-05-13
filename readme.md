[![Bioconductor Time](http://bioconductor.org/shields/years-in-bioc/categoryCompare.svg)](http://bioconductor.org/packages/release/bioc/html/categoryCompare.html "Bioconductor status")

[![Bioconductor Availability](http://bioconductor.org/shields/availability/release/categoryCompare.svg)](http://bioconductor.org/packages/release/bioc/html/categoryCompare.html#archives "Platform availability") 
[![Bioconductor Downloads](http://bioconductor.org/shields/downloads/categoryCompare.svg)](http://bioconductor.org/packages/stats/bioc/categoryCompare.html "Percentile downloads")
[![Bioconductor Commits](http://bioconductor.org/shields/commits/bioc/categoryCompare.svg)](http://bioconductor.org/packages/release/bioc/html/categoryCompare.html#svn_source "svn commits")
[![Support posts](http://bioconductor.org/shields/posts/categoryCompare.svg)](https://support.bioconductor.org/t/categorycompare/ "Bioconductor support posts")

[![Build Status](https://travis-ci.org/rmflight/categoryCompare.svg?branch=master)](https://travis-ci.org/rmflight/categoryCompare "travis build status") [![Bioconductor Release Build](http://bioconductor.org/shields/build/release/bioc/categoryCompare.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/categoryCompare/ "Bioconductor release build") [![Bioconductor Devel Build](http://bioconductor.org/shields/build/devel/bioc/categoryCompare.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/categoryCompare/ "Bioconductor devel build")

# categoryCompare

A Bioconductor package for meta analysis of high-throughput datasets using 
enriched feature annotations instead of just the features themselves.

See the "Description" file for additional requirements.

## Documentation

The [Vignette][vignLink] provides a description of the thinking behind
this package as well as a toy example for demonstration purposes.

## Installation

Installation of this package from Github requires the [devtools][devtoolsLink]
package.

```r
install.packages("devtools")
library(devtools)
install_github("categoryCompare", "rmflight")
```

Alternatively, you can install **categoryCompare** from Bioconductor itself.

```r
source("http://bioconductor.org/biocLite.R")
biocLite("categoryCompare")
```


[vignLink]: http://rmflight.github.io/categoryCompare/index.html "categoryCompare Vignette"
[devtoolsLink]: https://github.com/hadley/devtools "devtools"

## Citation

Flight RM, Harrison BJ, Mohammad F, Bunge MB, Moon LDF, Petruska JC and Rouchka EC (2014). .CATEGORYCOMPARE, an analytical tool based on feature annotations.
_Frontiers in Genetics_. [link](http://dx.doi.org/10.3389/fgene.2014.00098)
