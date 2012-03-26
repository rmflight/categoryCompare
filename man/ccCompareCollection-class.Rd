\name{ccCompareCollection-class}
\Rdversion{1.1}
\docType{class}
\alias{ccCompareCollection-class}
\alias{ccCompareCollection}

\title{Class \code{"ccCompareCollection"}}
\description{
Holds multiple \code{ccCompareResult} objects.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ccCompareCollection", ...)}.

These are not normally created by the user, but rather by \code{\link{ccCompare}} while performing the categorical comparisons for each type of category
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"}  }
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
No methods defined with class "ccCompareCollection" in the signature.
}

\author{
Robert M Flight
}

\seealso{
\code{\linkS4class{ccCompareResult}} \code{\link{ccCompare}}
}
\examples{
showClass("ccCompareCollection")
}
\keyword{classes}
