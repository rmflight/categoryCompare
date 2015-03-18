\name{HyperGParamsCC-class}
\Rdversion{1.1}
\docType{class}
\alias{HyperGParamsCC-class}
\alias{HyperGParamsCC}
\alias{GOHyperGParams}
\alias{GOHyperGParams-class}
\alias{GOHyperGParamsCC-class}
\alias{GOHyperGParamsCC}
\alias{KEGGHyperGParamsCC-class}
\alias{KEGGHyperGParamsCC}
\alias{KEGGHyperGParams-class}
\alias{KEGGHyperGParams}
\alias{fdr,HyperGParamsCC-method}
\alias{fdr<-,HyperGParamsCC-method}

\title{Class "HyperGParamsCC"}
\description{
This class extends the \code{HyperGParams} class in \code{Category} by providing options for multiple testing and the storing of extra data in addition to the gene list of interest (not currently used, but might be in the future).
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("HyperGParamsCC", ...)}. In general the user will not create these directly, but they are created and used by \code{} to carry out the enrichment calculations.
}
\section{Slots}{
  \describe{
    \item{\code{fdr}:}{Object of class \code{"numeric"} The number of FDR runs to perform }
    \item{\code{data}:}{Object of class \code{"data.frame"} Extra data stored in the object }
    \item{\code{geneIds}:}{Object of class \code{"ANY"} The genes of interest }
    \item{\code{universeGeneIds}:}{Object of class \code{"ANY"} The gene universe or background used (all the genes on the chip) }
    \item{\code{annotation}:}{Object of class \code{"character"} The annotation package used to get information about the geneIds }
    \item{\code{datPkg}:}{Object of class \code{"DatPkg"} Generated automatically from the \code{annotation} slot }
    \item{\code{categorySubsetIds}:}{Object of class \code{"ANY"} A specific set of category IDs that one wants to restrict the testing to }
    \item{\code{categoryName}:}{Object of class \code{"character"} What type of category to use, currently either "GO" or "KEGG" }
    \item{\code{pvalueCutoff}:}{Object of class \code{"numeric"} What should be the p-value to decide significance }
    \item{\code{testDirection}:}{Object of class \code{"character"} "over" or "under" represented annotation terms }
  }
}
\section{Extends}{
Class \code{"\linkS4class{GOHyperGParams}"}, directly.
}
\section{Methods}{
No methods defined with class "HyperGParamsCC" in the signature.
}
\author{
Robert M Flight
}
\seealso{
%%  objects to See Also as \code{\link{fun}}, ~
\code{\linkS4class{HyperGResultCC}} \code{\link{ccEnrich}} \code{Category-package}
}
\examples{
showClass("HyperGParamsCC")
}
\keyword{classes}
