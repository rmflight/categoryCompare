\name{ccEnrich-method}
\alias{ccEnrich-method}
\alias{ccEnrich}
\alias{ccEnrich,ccGeneList-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform annotation enrichment for multiple gene lists
}
\description{
Takes a \code{\link{ccGeneList}} object containing all the information needed to perform enrichment calculations for Gene Ontology.
}
\usage{
ccEnrich(ccGeneList)
}

\arguments{
  \item{ccGeneList}{
  A \code{\link{ccGeneList}} object, which is really just a list of lists, with some extra slots to tell us how to examine results. Each entry in the list should be named to allow identification later on. Each sub list should contain a vector \env{genes} denoting the genes of interest, a vector \env{universe} denoting the gene background (i.e. all genes on the chip), and an entry \env{annotation} denoting an organism database package (such as \pkg{org.Hs.eg.db}). See \code{\linkS4class{ccGeneList}} for more details regarding this object.
}
}
\details{
  This function is essentially a wrapper for \code{\link{hyperGTestCC}} that performs all of the calculations for the many gene lists in one go, returning a list of \code{\link{HyperGResultCC}} objects, one for each of the \code{ccTypes} and each gene list. These various \code{\link{HyperGResultCC}} objects can then be accessed and results compared among the lists for each of the ontologies
  }
\value{
A list of \code{\link{HyperGResultCC}} objects, one for each \code{ccType} and gene list, returned as \code{\link{ccEnrichResult}} objects for each \code{ccType}. This can be passed with a \code{\link{ccOptions}} object to \code{ccCompare} to generate actual annotation comparisons.
}

\author{
Robert M Flight
}

\seealso{
\code{\link{ccGeneList}}, \code{\link{hyperGTestCC}}, \code{\link{ccEnrichResult}}
}
\examples{
\dontrun{
require(GO.db)
require(KEGG.db)
require(org.Hs.eg.db)
}
data(ccData)

g10 <- (unique(table10$Entrez[1:100]))
g48 <- (unique(table48$Entrez[1:100]))

list10 <- list(genes=g10, universe=gUniverse, annotation="org.Hs.eg.db")
list48 <- list(genes=g48, universe=gUniverse, annotation="org.Hs.eg.db")

geneLists <- list(T10=list10, T48=list48)
geneLists <- new("ccGeneList", geneLists, ccType=c("BP","KEGG"))
geneLists <- new("ccGeneList", geneLists, ccType=c("CC","KEGG"))

# set number of fdr runs to 0 to speed up runtime, not generally recommended.
geneLists <- new("ccGeneList", geneLists, ccType = c('BP','KEGG'), pvalueCutoff=0.01, fdr=0)
# enrichLists <- ccEnrich(geneLists)
}
