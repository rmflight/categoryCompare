\name{ccCompareResult-class}
\Rdversion{1.1}
\docType{class}
\alias{ccCompareResult-class}
\alias{ccCompareResult}
\alias{mainGraph,ccCompareResult-method}
\alias{mainTable,ccCompareResult-method}
\alias{allAnnotation,ccCompareResult-method}
\alias{mainGraph}
\alias{mainTable}
\alias{allAnnotation}

\title{Class \code{"ccCompareResult"}}
\description{
Holds the results from a single category comparison
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ccCompareResult", ...)}.
%%   describe objects here  
}
\section{Slots}{
  \describe{
    \item{\code{mainGraph}:}{Object of class \code{"graph"}. Holds the graph describing the relationships between the annotations }
    \item{\code{subGraph}:}{Object of class \code{"list"}. Not currently used }
    \item{\code{mainTable}:}{Object of class \code{"data.frame"}. Table of results, with all the various statistics for each annotation in the category }
    \item{\code{allAnnotation}:}{Object of class \code{"list"}. For each annotation, which genes from which comparison are annotated to that particular annotation }
    \item{\code{categoryName}:}{Object of class \code{"character"}. Which category (e.g. GO, KEGG, etc) was used }
    \item{\code{ontology}:}{Object of class \code{"character"}. If GO, which ontology was used }
  }
}
\section{Methods}{
  \describe{
    \item{allAnnotation}{\code{signature(object = "ccCompareResult")}: ... }
    \item{mainGraph}{\code{signature(object = "ccCompareResult")}: ... }
    \item{mainTable<-}{\code{signature(object = "ccCompareResult")}: ... }
	 }
}

\author{
Robert M Flight
}

\seealso{
\code{\link{ccCompare}} \code{\linkS4class{ccCompareCollection}}
}
\examples{
showClass("ccCompareResult")
}
\keyword{classes}
