\name{ccOutCyt-methods}
\docType{methods}
\alias{ccOutCyt-methods}
\alias{ccOutCyt}
\alias{ccOutCyt,ccCompareResult,ccOptions-method}
\title{  Methods for Function \code{ccOutCyt} in Package \pkg{categoryCompare} }
\description{
Passes a \code{ccCompareResult} object to \code{Cytoscape} for interactive visualization of \code{ccCompare} results.
}

\section{Methods}{
\describe{

\item{\code{signature(ccCompRes = "ccCompareResult", ccOpts = "ccOptions", ...)}}{
At a minimum, this method requires a \code{ccCompareResult} and a \code{ccOptions} to work.
}

... may include:

\item{layout = "character"}{to override the default layout set by \code{ccCompare}, as well as options}
\item{postText = "character"}{to add a user set string to the Cytoscape window}

In addition, any of the arguments to \code{CytoscapeWindow} may also be set, such as \code{host} or \code{port}.
}
}
\details{
Note that only some basic, required methods have been imported from \code{RCy3} for use with \code{categoryCompare}, and these are hidden in the functions within \code{categoryCompare} and are not visible to the user. If access to all the functionality of \code{RCytoscape} is desired (and trust me, there is a lot of useful stuff in there), then the user should use \code{library(RCy3)} directly.

It should also be noted that deletion of edges via \code{RCy3} is slow, so some edge filtering should be done by \code{\link{breakEdges}} prior to using \code{ccOutCyt}.
}
\seealso{
\code{\linkS4class{ccCompareResult}} \code{\linkS4class{ccOptions}} \code{\link{ccCompare}} \code{CytoscapeWindowClass}
}
\examples{
\dontrun{
hasCy <- (if (.Platform$OS.type \%in\% "windows") { (length(grep("Cytoscape", system("tasklist", intern=TRUE))) > 0)})

if hasCy {

	ccResults$BP <- breakEdges(ccResults$BP, 0.8)
	cwObj <- ccOutCyt(ccResults$BP,ccOpts)
	Sys.sleep(10)
	RCytoscape::deleteWindow(cwObj)
	}
}}
\keyword{methods}
