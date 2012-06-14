\name{breakEdges-methods}
\docType{methods}
\alias{breakEdges-methods}
\alias{breakEdges,ccCompareResult,numeric-method}
\alias{breakEdges,CytoscapeWindowClass,numeric-method}
\title{  Methods for Function \code{breakEdges} in Package \pkg{categoryCompare} }
\description{
Methods for function \code{breakEdges} in package \pkg{categoryCompare}
}
\section{Methods}{
\describe{

\item{\code{signature(cwObject = "ccCompareResult", cutoff = "numeric")}}{
Allows one to remove edges in the \code{ccCompareResult} \code{mainGraph} slot prior to passing it into Cytoscape for visualization. Given that the number of edges can be rather large (especially for Gene Ontology) this can easily speed up the transfer, without actually losing any information.
}

\item{\code{signature(cwObject = "CytoscapeWindowClass", cutoff = "numeric")}}{
Once an annotation graph is in Cytoscape, remove edges above or below the cutoff. Note that this does not affect the original graph in the \code{ccCompareResult} object.
}
}}

\author{
Robert M Flight
}

\seealso{
\code{\linkS4class{CytoscapeWindowClass}} \code{\link{breakEdges}} \code{\linkS4class{ccCompareResult}} \code{\link{ccOutCyt}}
}
\examples{
data(ccData)

# breaking the edges in a ccCompareResult
ccResults$CC <- breakEdges(ccResults$CC, 0.8)
\dontrun{
hasCy <- (if (.Platform$OS.type \%in\% "windows") { (length(grep("Cytoscape", system("tasklist", intern=TRUE))) > 0)})

if hasCy {
	
	cwObj <- ccOutCyt(ccResults$CC,ccOpts)
	# now breaking them in the CytoscapeWindow object
	breakEdges(cwObj, 0.85)
	Sys.sleep(10)
	RCytoscape::deleteWindow(cwObj)
	}
}}
\keyword{methods}
\keyword{  other possible keyword(s)  }
