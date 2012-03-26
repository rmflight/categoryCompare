\name{ccGeneList-class}
\Rdversion{1.1}
\docType{class}
\alias{ccGeneList-class}
\alias{ccGeneList}
\alias{listNames,ccGeneList-method}
\alias{fdr,ccGeneList-method}
\alias{fdr<-,ccGeneList-method}
\alias{pvalueCutoff,ccGeneList-method}
\alias{pvalueCutoff<-,ccGeneList-method}
\alias{ccType,ccGeneList-method}
\alias{ccType<-,ccGeneList-method}
\alias{testDirection,ccGeneList-method}
\alias{ccType}
\alias{ccType<-}

\title{Class "ccGeneList"}
\description{
This stores the actual gene lists and related information that will be used in categoryCompare.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ccGeneList", list))}. \code{ccGeneList} is actually just an extension of R list objects. The input \code{list} should be a list of lists. See \code{Details} for more information.}
\section{Slots}{
  \describe{
    \item{\code{fdr}:}{Object of class \code{"numeric"} The number of \code{fdr} runs to perform to account for different list sizes and term dependence }
    \item{\code{pvalueCutoff}:}{Object of class \code{"numeric"} Value used to determine whether or not a particular term is significant or not }
    \item{\code{ccType}:}{Object of class \code{"character"} What types of annotations to use. Currently supported ones include "BP", "MF", "CC" (from Gene Ontology) and "KEGG" }
    \item{\code{testDirection}:}{Object of class \code{"character"} Are you interested in "over" or "under" represented annotations }
  }
}
\section{Methods}{
  \describe{
    \item{fdr}{\code{signature(object = "ccGeneList")}: how many random runs to perform }
    \item{fdr<-}{\code{signature(object = "ccGeneList")}: change the number of random runs }
    \item{pvalueCutoff}{\code{signature(object = "ccGeneList")}: what is the pvalue to consider significant }
    \item{pvalueCutoff<-}{\code{signature(object = "ccGeneList")}: change the cutoff for significance }
    \item{ccType}{\code{signature(object = "ccGeneList")}: what type of annotations are going to be examined }
    \item{ccType<-}{\code{signature(object = "ccGeneList")}: change the type of annotations to examine }
    \item{testDirection}{\code{signature(object = "ccGeneList")}: query for "over" or "under" represented annotations }
    \item{testDirection<-}{\code{signature(object = "ccGeneList")}: change the type of representation ("over" or "under") }
    \item{listNames}{\code{signature(object = "ccGeneList")}: what are the names of the lists contained }
   }
}
\section{Details}{
  The input list should be a list of lists, with at least three sub-lists. 
  
  \code{testList <- list(list1=list(genes='...',universe='...',annotation='...'), list2=list(...))}
  \describe{
  \item{genes}{: These are the gene identifiers of the genes that are of interest (differentially expressed genes)}
  \item{universe}{: All of the genes that were measured in this particular experiments (i.e. all the genes on the chip)}
  \item{annotation}{: What organism or chip do these ID's come from (e.g. "org.Hs.eg.db" for Human Entrez gene ID's, "hgu133a.db" for probe ID's from the Affymetrix U133A chip)}
  \item{data}{: A data-frame that contains extra information about the genes of interest. At the very least, the data-frame must have a column \code{ID} that matches the ID's contained in \code{genes}}
  }
  
  What actually happens when running \code{ccEnrich} is that the appropriate \code{HyperGParamsCC} objects are generated for each geneList and each type of annotation (e.g. BP, CC, KEGG), and then the calculations performed on each one.
}

\author{
Robert M Flight
}
\note{
The \code{ccGeneList} object is what will undergo all of the enrichment calculations. When the results are combined with the \code{\link{ccOptions}} object, we can get our results of actual comparisons between experiments.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\linkS4class{ccOptions}}
}
\examples{
data(ccData)
g10 <- (unique(table10$Entrez[1:100]))
g48 <- (unique(table48$Entrez[1:100]))

list10 <- list(genes=g10, universe=gUniverse, annotation="org.Hs.eg.db")
list48 <- list(genes=g48, universe=gUniverse, annotation="org.Hs.eg.db")

geneLists <- list(T10=list10, T48=list48)
geneLists <- new("ccGeneList", geneLists, ccType=c("BP","KEGG"))
geneLists
}
\keyword{classes}