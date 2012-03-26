\name{minCount}
\alias{minCount}
\alias{minCount<-}
\title{
minCount
}
\description{
Extracts and sets the minimum number of genes that an annotation must have to be considered in subsequent steps.
}
\usage{
	minCount(object)
}
\arguments{
	\item{object}{ This will be either a \code{HyperGResultCC}, \code{ccEnrichResult}, or \code{ccEnrichCollection} object. See Details for more information.}
}
\details{
	\code{minCount(object)} fetches the set \code{minCount} for \code{HyperGResultCC} and \code{ccEnrichResult} objects
	
	\code{minCount(object)<-} will set the \code{minCount} for \code{HyperGResultCC} objects, and when applied to \code{ccEnrichResult} and \code{ccEnrichCollection} sets the \code{minCount} for all of the contained objects, so be careful if you want to use different \code{minCounts} for different results
}
\author{
Robert M Flight
}
\seealso{
\code{\linkS4class{HyperGResultCC}} \code{\linkS4class{ccEnrichResult}} \code{\linkS4class{ccEnrichCollection}}
}
\examples{
data(ccData)
enrichLists
minCount(enrichLists) <- 5
enrichLists
}
