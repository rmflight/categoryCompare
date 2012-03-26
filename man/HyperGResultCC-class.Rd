\name{HyperGResultCC-class}
\Rdversion{1.1}
\docType{class}
\alias{HyperGResultCC}
\alias{HyperGResultCC-class}
\alias{GOHyperGResultCC-class}
\alias{KEGGHyperGResultCC-class}
\alias{GOHyperGResultCC}
\alias{KEGGHyperGResultCC}
\alias{fdr,HyperGResultCC-method}
\alias{fdrvalues,HyperGResultCC-method}
\alias{pCC,HyperGResultCC-method}
\alias{pvalueCutoff<-,HyperGResultCC-method}
\alias{pvalueType,HyperGResultCC-method}
\alias{pvalueType<-,HyperGResultCC-method}
\alias{fdrvalues-method}
\alias{fdrvalues}
\alias{minCount,HyperGResultCC-method}
\alias{minCount<-,HyperGResultCC-method}

\title{Class "HyperGResultCC"}
\description{
Contains the results of performing a hypergeometric test on a \code{\linkS4class{HyperGParams}} object. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("HyperGResultCC", ...)}.
%%   describe objects here  
}
\section{Slots}{
  \describe{
    \item{\code{fdr}:}{Object of class \code{"numeric"} The number of FDR runs performed }
    \item{\code{fdrvalues}:}{Object of class \code{"numeric"} The FDR values generated }
    \item{\code{pvalueType}:}{Object of class \code{"character"} Whether to use p-values or FDR values in determing the significant terms returned }
    \item{\code{data}:}{Object of class \code{"data.frame"} Extra data }
    \item{\code{pvalues}:}{Object of class \code{"numeric"} P-values calculated for each term }
    \item{\code{oddsRatios}:}{Object of class \code{"numeric"}  }
    \item{\code{expectedCounts}:}{Object of class \code{"numeric"}  }
    \item{\code{catToGeneId}:}{Object of class \code{"list"}  }
    \item{\code{organism}:}{Object of class \code{"character"}  }
    \item{\code{annotation}:}{Object of class \code{"character"}  }
    \item{\code{geneIds}:}{Object of class \code{"ANY"}  }
    \item{\code{testName}:}{Object of class \code{"character"}  }
    \item{\code{pvalueCutoff}:}{Object of class \code{"numeric"}  }
    \item{\code{testDirection}:}{Object of class \code{"character"}  }
  }
}
\section{Extends}{
Class \code{"\linkS4class{HyperGResult}"}, directly.
Class \code{"\linkS4class{HyperGResultBase}"}, by class "HyperGResult", distance 2.
}
\section{Methods}{
  \describe{
    \item{fdr}{\code{signature(object = "HyperGResultCC")}: ... }
    \item{fdrvalues}{\code{signature(object = "HyperGResultCC")}: ... }
    \item{pCC}{\code{signature(object = "HyperGResultCC")}: ... }
    \item{pvalueCutoff<-}{\code{signature(r = "HyperGResultCC")}: ... }
    \item{pvalueType}{\code{signature(object = "HyperGResultCC")}: ... }
    \item{pvalueType<-}{\code{signature(object = "HyperGResultCC")}: ... }
    \item{minCount}{\code{signature(object = "HyperGResultCC")}: ...}
    \item{minCount<-}{\code{signature(object = "HyperGResultCC")}: ...}
	 }
}

\author{
Robert M Flight
}
\seealso{
	\code{\link{hyperGTestCC}}
}
\examples{
showClass("HyperGResultCC")
}
\keyword{classes}
