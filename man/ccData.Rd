\name{ccData}
\alias{ccData}
\alias{gseaRes}
\alias{gUniverse}
\alias{table10}
\alias{table48}
\alias{enrichLists}
\alias{ccResults}
\alias{ccResultsBPHier}
\alias{ccOpts}
\alias{geneLists}
\docType{data}
\title{
Test data for \code{categoryCompare}
}
\description{
Processed data from the \code{estrogen} example data set
}
\usage{data(ccData)}
\format{
  \code{table10}{: Log-ratio output from \pkg{limma} for the comparison of presence-absence of estrogen at 10 hours}
  
  \code{table48}{: Log-ratio output from \pkg{limma} for the comparison of presence-absence of estrogen at 48 hours}
  
  \code{gUniverse}{: All of the genes measured on the chip}

	\code{gseaRes}{: Toy results of GSEA analysis of 3 different tissues}
	
	\code{enrichLists}{: Apply \code{\link{ccEnrich}} to a ccGeneList from \code{table10} and \code{table48}}
	
	\code{ccResults}{: Apply \code{\link{ccCompare}} to \code{enrichLists}}
	
	\code{ccResultsHier}{: Modify \code{enrichLists$CC} to use a "hierarchical" layout}
	
	\code{geneLists}{: a \code{ccGeneList} generated from genes in table10 and table48}
	
	\code{ccOpts}{: a \code{ccOptions} object describing what we are going to do as far as feature list comparisons}
	
}
\source{
Taken from the \pkg{estrogen} package in Bioconductor, and then processed using the normal \pkg{affy} and \pkg{limma} tools.
}
\author{
Robert M Flight
}
\seealso{
\code{\linkS4class{ccGeneList}} \code{\linkS4class{ccEnrichCollection}} \code{\linkS4class{ccCompareCollection}} \code{\link{ccEnrich}} \code{\link{ccCompare}}
}
\examples{
data(ccData)

}

\keyword{datasets}
