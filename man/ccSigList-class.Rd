\name{ccSigList-class}
\Rdversion{1.1}
\docType{class}
\alias{ccSigList}
\alias{ccSigList-class}
\alias{annotation,ccSigList-method}
\alias{category,ccSigList-method}
\alias{ontology,ccSigList-method}
\alias{sigID,ccSigList-method}

\title{Class \code{"ccSigList"}}
\description{
Holds a generic list of significant annotations. Allows one to use Bioconductor annotation packages, or when combined into a \code{GENccEnrichResult}, use custom annotation / gene mappings.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ccSigList", ...)}.
%%   describe objects here  
}
\section{Slots}{
  \describe{
    \item{\code{sigID}:}{Object of class \code{"character"}  }
    \item{\code{categoryName}:}{Object of class \code{"character"}  }
    \item{\code{ontology}:}{Object of class \code{"character"}  }
    \item{\code{annotation}:}{Object of class \code{"character"}  }
  }
}
\section{Methods}{
  \describe{
    \item{annotation}{\code{signature(object = "ccSigList")}: ... }
    \item{category}{\code{signature(object = "ccSigList")}: ... }
    \item{ontology}{\code{signature(object = "ccSigList")}: ... }
    \item{sigID}{\code{signature(object = "ccSigList")}: ... }
	 }
}

\author{
Robert M Flight
}

\seealso{
\code{\linkS4class{GENccEnrichResult}} \code{\link{ccCompareGeneric}} 
}
\examples{
showClass("ccSigList")
}
\keyword{classes}
