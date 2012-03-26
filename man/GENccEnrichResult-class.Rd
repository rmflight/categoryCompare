\name{GENccEnrichResult-class}
\Rdversion{1.1}
\docType{class}
\alias{GENccEnrichResult-class}
\alias{GENccEnrichResult}
\alias{[,GENccEnrichResult,ANY,ANY-method}
\alias{category,GENccEnrichResult-method}
\alias{ccCompareGeneric,GENccEnrichResult,ccOptions-method}
\alias{geneAnnMapping,GENccEnrichResult-method}
\alias{graphType,GENccEnrichResult-method}
\alias{graphType<-,GENccEnrichResult-method}
\alias{ontology,GENccEnrichResult-method}
\alias{geneAnnMapping}
\alias{sigID}

\title{Class \code{"GENccEnrichResult"}}
\description{
Holds generic \code{ccEnrich} type results
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GENccEnrichResult", ...)}.
%%   describe objects here  
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} The actual list containing the \code{ccEnrichResults} }
    \item{\code{categoryName}:}{Object of class \code{"character"}  }
    \item{\code{ontology}:}{Object of class \code{"character"}  }
    \item{\code{geneAnnMapping}:}{Object of class \code{"namedList"}  }
    \item{\code{graphType}:}{Object of class \code{"character"}  }
    \item{\code{names}:}{Object of class \code{"character"}  }
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
    \item{[}{\code{signature(x = "GENccEnrichResult", i = "ANY", j = "ANY")}: Subsets the object to just those lists that are desired }
    \item{categoryName}{\code{signature(object = "GENccEnrichResult")}:  }
    \item{ccCompareGeneric}{\code{signature(gccResult = "GENccEnrichResult", ccOptions = "ccOptions")}: ... }
    \item{geneAnnMapping}{\code{signature(object = "GENccEnrichResult")}: ... }
    \item{graphType}{\code{signature(object = "GENccEnrichResult")}: ... }
    \item{graphType<-}{\code{signature(object = "GENccEnrichResult")}: ... }
    \item{ontology}{\code{signature(object = "GENccEnrichResult")}: ... }
	 }
}

\author{
Robert M Flight
}

\seealso{
\code{\link{ccCompareGeneric}} \code{\link{ccSigList}}
}
\examples{
data(ccData)
locA <- grep("A",gseaRes$Tissues)
locL <- grep("L",gseaRes$Tissues)
locM <- grep("M",gseaRes$Tissues)

A <- new("ccSigList",sigID=gseaRes$KEGGID[locA],categoryName="KEGG",annotation="org.Mm.eg")
L <- new("ccSigList",sigID=gseaRes$KEGGID[locL],categoryName="KEGG",annotation="org.Mm.eg")
M <- new("ccSigList",sigID=gseaRes$KEGGID[locM],categoryName="KEGG",annotation="org.Mm.eg")
ccEnrichCol <- list(A=A,L=L,M=M)
ccEnrichCol <- new("GENccEnrichResult", ccEnrichCol, categoryName="KEGG")
}
\keyword{classes}
