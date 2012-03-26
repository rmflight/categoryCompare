\name{mergedData-class}
\Rdversion{1.1}
\docType{class}
\alias{mergedData-class}

\title{Class \code{"mergedData"}}
\description{
Stores merged data tables from the "data" entry in a \code{ccGeneList}. This is useful for output later.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("mergedData", ...)}.
%%   describe objects here  
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"}  }
    \item{\code{useIDName}:}{Object of class \code{"character"}  }
    \item{\code{names}:}{Object of class \code{"character"}  }
    \item{\code{row.names}:}{Object of class \code{"data.frameRowLabels"}  }
    \item{\code{.S3Class}:}{Object of class \code{"character"}  }
  }
}
\section{Extends}{
Class \code{"\linkS4class{data.frame}"}, directly.
Class \code{"\linkS4class{list}"}, by class "data.frame", distance 2.
Class \code{"\linkS4class{oldClass}"}, by class "data.frame", distance 2.
Class \code{"\linkS4class{data.frameOrNULL}"}, by class "data.frame", distance 2.
Class \code{"\linkS4class{vector}"}, by class "data.frame", distance 3.
}
\section{Methods}{
\code{signature(saveObj = "list", compareResult = "ccCompareResult", mergedData = "mergedData")}
}

\author{
Robert M. Flight
}

\seealso{
\code{\link{mergeLists}} \code{\link{cytOutData}}
}
\examples{
showClass("mergedData")

data(ccData)
mergeDat <- mergeLists(geneLists, ccOpts)

}
\keyword{classes}
