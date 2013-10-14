\name{ccEnrichResult-class}
\Rdversion{1.1}
\docType{class}
\alias{ccEnrichResult}
\alias{ccEnrichResult-class}
\alias{GOccEnrichResult-class}
\alias{GOccEnrichResult}
\alias{KEGGccEnrichResult-class}
\alias{KEGGccEnrichResult}
\alias{fdr,ccEnrichResult-method}
\alias{pvalueCutoff,ccEnrichResult-method}
\alias{pvalueCutoff<-,ccEnrichResult-method}
\alias{pvalueType<-,ccEnrichResult-method}
\alias{pvalueType,ccEnrichResult-method}
\alias{minCount,ccEnrichResult-method}
\alias{minCount<-,ccEnrichResult-method}
\alias{graphType,ccEnrichResult-method}
\alias{graphType<-,ccEnrichResult-method}
\alias{[,ccEnrichResult,ANY,ANY-method}
\alias{[,ccEnrichResult,ANY,ANY,ANY-method}

\title{Class "ccEnrichResult"}
\description{
Acts as a container object for multiple \code{\linkS4class{HyperGResultCC}} objects.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ccEnrichResult", ...)}.
%%   describe objects here  
}

\section{Extends}{
Class \code{"\linkS4class{namedList}"}, directly.
Class \code{"\linkS4class{list}"}, by class "namedList", distance 2.
Class \code{"\linkS4class{vector}"}, by class "namedList", distance 3.
Class \code{"\linkS4class{AssayData}"}, by class "namedList", distance 3.
}
\section{Methods}{
  \describe{
    \item{fdr}{\code{signature(object = "ccEnrichResult")}: get the number of runs using random feature lists were performed }
    \item{pvalueCutoff}{\code{signature(r = "ccEnrichResult")}: what is the pvalueCutoff to determine significant annotations }
    \item{pvalueCutoff<-}{\code{signature(r = "ccEnrichResult")}: change the pvalueCutoff for an annotation to be considered significant }
    \item{pvalueType<-}{\code{signature(object = "ccEnrichResult")}: change whether p-values used are from "FDR" or raw p-values }
    \item{minCount}{\code{signature(object = "ccEnrichResult")}: how many features need to belong to an annotation to be reported }
    \item{minCount<-}{\code{signature(object = "ccEnrichResult")}: adjust the minCount}
    \item{graphType}{\code{signature(object = "ccEnrichResult")}: what type of graph should be generated (generally set by the class of object)}
    \item{graphType<-}{\code{signature(object = "ccEnrichResult")}: change the type of graph to generate by \code{ccCompare}}
	 }
}
\author{
Robert M Flight
}
\examples{
data(ccData)
enrichRes <- enrichLists[[1]]
fdr(enrichRes)
pvalueType(enrichRes)
enrichRes
pvalueType(enrichRes) <- 'pval'
enrichRes

pvalueCutoff(enrichRes)
pvalueCutoff(enrichRes) <- 0.01
enrichRes
}
\keyword{classes}
