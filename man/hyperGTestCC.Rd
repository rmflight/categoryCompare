\name{hyperGTestCC}
\alias{hyperGTestCC}
\alias{hyperGTestCC,HyperGParamsCC-method}
\alias{hyperGTestCC,HyperGParamsCC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypergeometric testing with false discovery rate
}
\description{
Performs the hypergeometric testing for \code{\link{HyperGParamsCC}} objects.
}
\usage{
hyperGTestCC(p)
}

\arguments{
  \item{p}{ A \code{\link{HyperGParamsCC}} object
}
}

\details{
	This is the heart of categoryCompare, the function that calculates the HyperGeometric statistics for the given categories of annotation for each gene list.
}

\value{ Returns a \code{\link{HyperGResultCC}} object
}
\author{
Robert M Flight
}
\seealso{
\code{\link{HyperGParamsCC}} \code{\link{HyperGResultCC}} \code{\link{GOHyperGParamsCC}} \code{\link{KEGGHyperGParamsCC}} \code{\link{GOHyperGResultCC}} \code{\link{KEGGHyperGResultCC}}
}
\examples{
require(GO.db)
require(org.Hs.eg.db)
data(ccData)
g10 <- unique(table10$Entrez)
testGO <- new("GOHyperGParamsCC", geneIds=g10, universeGeneIds=gUniverse, 
annotation="org.Hs.eg.db", ontology="CC", conditional=FALSE, 
testDirection="over",fdr=0, pvalueCutoff = 0.01)
# ccHypRes <- hyperGTestCC(testGO)
# summary(ccHypRes)
}
