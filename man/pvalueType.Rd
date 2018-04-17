\name{pvalueType}
\alias{pvalueType}
\alias{pvalueType<-}
\title{
Type of p-values to return from object
}
\description{
Queries or sets the type of p-values to return from objects, either base calculated (pvals) or from fdr calculations (fdr)
}
\usage{
	pvalueType(object)
}
\arguments{
	\item{object}{ Can be \code{HyperGResultCC}, \code{ccEnrichResult}, \code{ccEnrichCollection}. See Details for more information}
}
\details{
	\code{pvalueType(object)} gets the type of \code{p-values} to be returned from \code{HyperGResultCC} and \code{ccEnrichResult} objects

	\code{pvalueType(object)<-} will set the type of \code{p-values} to be returned from \code{HyperGResultCC}, \code{ccEnrichResult}, \code{ccEnrichCollection}. Note that for a \code{ccEnrichCollection}, the type is changed for all contained \code{ccEnrichResults}
}
\author{
Robert M Flight
}
\seealso{
\code{\linkS4class{HyperGResultCC}} \code{\linkS4class{ccEnrichResult}} \code{\linkS4class{ccEnrichCollection}}
}

\examples{
# pvalueType-Methods
data(ccData)

\dontrun{pvalueType(enrichLists) # this returns an error}
pvalueType(enrichLists[[1]])
pvalueType(enrichLists[[1]][[1]])

# change the type for one of the results
pvalueType(enrichLists[[1]]) <- 'pval'  # Not recommended practice
enrichLists

# change for all of the results
pvalueType(enrichLists) <- 'pval'
enrichLists
}
