\name{fdr}
\alias{fdr}
\alias{fdr<-}
\title{
Number of FDR runs to perform
}
\description{
Queries or sets the number of random runs to perform to generate an estimate of the false discovery rate. Defaults to 50
}
\usage{
	fdr(object)
}
\arguments{
	\item{object}{ Can be \code{ccGeneList}, \code{HyperGParamsCC}, \code{HyperGResultCC}, \code{ccEnrichResult} See Details for more information.}
}
\details{
	\code{fdr(object)} gets the number of \code{fdr} runs for \code{ccGeneList}, \code{HyperGParamsCC}, \code{HyperGResultCC}, \code{ccEnrichResult}
	
	\code{fdr(object)<-} will set the number of \code{fdr} runs to be used by \code{ccEnrich} and \code{HyperGTestCC} when performing calculations on either a \code{ccGeneList} or \code{HyperGParamsCC}, respectively
}
\author{
Robert M Flight
}
\seealso{
\code{\linkS4class{HyperGResultCC}} \code{\linkS4class{ccEnrichResult}} \code{\link{ccGeneList}} \code{\link{HyperGParamsCC}}
}

