\name{ccEnrichCollection-class}
\Rdversion{1.1}
\docType{class}
\alias{ccEnrichCollection-class}
\alias{ccEnrichCollection}
\alias{pvalueCutoff<-,ccEnrichCollection-method}
\alias{pvalueType<-,ccEnrichCollection-method}
\alias{minCount<-,ccEnrichCollection-method}
\alias{graphType,ccEnrichCollection-method}
\alias{graphType<-,ccEnrichCollection-method}

\title{Class "ccEnrichCollection"}
\description{
Holds multiple classes of \code{\linkS4class{ccEnrichResult}} in one object to allow \code{\link{ccCompare}} to work on only the one object and generate all of the results of a comparison. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ccEnrichCollection", ...)}.
%%   describe objects here  
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"}  }
    \item{\code{names}:}{Object of class \code{"character"} The names (generally GO ontologies or KEGG, but can be changed) of each set of results }
  }
}
\section{Extends}{
Class \code{"\linkS4class{namedList}"}, directly.
Class \code{"\linkS4class{list}"}, by class "namedList", distance 2.
Class \code{"\linkS4class{vector}"}, by class "namedList", distance 3.
Class \code{"\linkS4class{AssayData}"}, by class "namedList", distance 3.
}
\section{Methods}{
  \describe{
    \item{pvalueCutoff<-}{\code{signature(r = "ccEnrichCollection")}: Changes the \code{pvalueCutoff} to be used to decide significant annotations for all of the contained \code{ccEnrichResult} objects }
    \item{pvalueType<-}{\code{signature(object = "ccEnrichCollection")}: Changes whether to use \code{p-values} or \code{fdr} values to determine those annotations that are significant in all of the contained \code{ccEnrichResult} objects }
    \item{minCount<-}{\code{signature(object = "ccEnrichCollection")}: how many features have to be annotated to a term to be reported as significant}
    \item{graphType}{\code{signature(object = "ccEnrichCollection")}: Gets the type of graph that should be output for this collection}
	 }
}

\author{
Robert M Flight
}
\seealso{
	\code{\link{ccEnrich}} \code{\link{hyperGTestCC}} \code{\link{ccCompare}} \code{\linkS4class{ccEnrichResult}}
}
\examples{
data(ccData)
enrichLists
}
\keyword{classes}
