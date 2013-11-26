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

[vignLink]: http://bioconductor.org/packages/devel/bioc/vignettes/categoryCompare/inst/doc/categoryCompare_vignette.pdf "categoryCompare Vignette"
[devtoolsLink]: https://github.com/hadley/devtools "devtools"
