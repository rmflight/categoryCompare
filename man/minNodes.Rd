\name{minNodes}
\alias{minNodes}
\alias{minNodes,CytoscapeWindowClass,numeric-method}
\title{
Delete nodes with less than a certain number of genes annotated
}
\description{
Deletes from the graph those annotations with less than a certain number of genes
}
\usage{
	minNodes(cwObj, cutoff)
}
\arguments{
	\item{cwObj}{ a \code{CytoscapeWindowClass} object returned from \code{ccOutCyt}}
	\item{cutoff}{ the minimum number of genes that an annotation must have}
}

\author{
Robert M Flight
}
\seealso{
\code{CytoscapeWindowClass} \code{\link{ccOutCyt}}
}
\examples{
\dontrun{
hasCy <- (if (.Platform$OS.type \%in\% "windows") { (length(grep("Cytoscape", system("tasklist", intern=TRUE))) > 0)})

if hasCy {
	data(ccData)
	ccResults$BP <- breakEdges(ccResults$BP, 0.8)
	cwObj <- ccOutCyt(ccResults$BP,ccOpts)
	
	minNodes(cwObj, 5)

}}}
