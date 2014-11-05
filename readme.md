[![Build Status](https://travis-ci.org/rmflight/categoryCompare.svg?branch=master)](https://travis-ci.org/rmflight/categoryCompare)

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
