\name{resetColors-methods}
\docType{methods}
\alias{resetColors-methods}
\alias{resetColors}
\alias{resetColors,numeric,ccOptions-method}
\title{ resetColors }
\description{
If the color of particular nodes have been modified from the original color scheme in \code{ccOptions}, this will reset them
}
\section{Methods}{
\describe{

\item{\code{signature(cwObj = "numeric", ccOpts = "ccOptions")}}{
What \code{CytoscapeWindow} to apply this to, and what \code{ccOptions} to use for the color scheme.
}

\subsection{Optional Arguments}{
Note that optional arguments include \code{node.attribute.name} (default is 'fillcolor') and \code{mode} (default is 'lookup')

}
}
}

\author{
Robert M Flight
}
\seealso{
\code{\linkS4class{ccOptions}} \code{\link{setNodeColorMapping}}
}
\keyword{methods}
\keyword{  other possible keyword(s)  }
