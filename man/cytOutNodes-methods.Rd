\name{cytOutNodes-methods}
\docType{methods}
\alias{cytOutNodes}
\alias{cytOutNodes-methods}
\alias{cytOutNodes,character,numeric,list-method}
\alias{cytOutNodes,character,numeric,missing-method}
\title{  Methods for Function \code{cytOutNodes}  }
\description{
Allows export of currently selected nodes in the Cytoscape window for data export
}
\section{Methods}{
\describe{

\item{\code{signature(descStr = "character", cwObj = "numeric", saveObj = "list")}}{
\code{descStr} is a string describing the nodes that are currently selected, \code{cwObj} is the \code{CytoscapeWindow} that the nodes are in, and then \code{saveObj} is a previously generated \code{cytOutNodes} list, and is optional.
}
}}
\examples{
\dontrun{
hasCy <- (if (.Platform$OS.type \%in\% "windows") {
(length(grep("Cytoscape", system("tasklist", intern=TRUE))) > 0)})

if hasCy {
	ccResults$BP <- breakEdges(ccResults$BP, 0.8)
	cwObj <- ccOutCyt(ccResults$BP,ccOpts)
	# user selects some nodes in Cytoscape
	RCy3::selectNodes(cwObj, c("GO:0007017", "GO:0000226", "GO:0007051", "GO:0007052"))
	savedNodes <- cytOutNodes("random1", cwObj) # save them
	# and selects some other nodes
	RCy3::selectNodes(cwObj,
	c("GO:0071103", "GO:0034728", "GO:0006323", "GO:0030261", "GO:0006334"),
	preserve.current.selection=FALSE)
	savedNodes <- cytOutNodes("random2", cwObj, savedNodes)

	}
}
}
\keyword{methods}
\keyword{  other possible keyword(s)  }
