\name{mergeLists-methods}
\docType{methods}
\alias{mergeLists-methods}
\alias{mergeLists}
\alias{mergeLists,ccGeneList,ccOptions-method}
\title{ Function \code{mergeLists} in Package \pkg{categoryCompare}}
\description{
Merges the gene lists or the data tables from a \code{ccGeneList} object, providing a single table with all the input data, that can then be queried later, using \code{cytTableOut}
}
\usage{
mergeLists(ccGeneList,ccOptions,isGene=TRUE)
}
\arguments{
	\item{ccGeneList}{a \code{ccGeneList} object}
	\item{ccOptions}{a \code{ccOptions} object}
	\item{isGene}{are the identifiers genes, or something else (metabolites, etc)}
}
\value{
	A \code{mergedData} object which is really just a glorified data frame. If the \code{ccGeneList} input had a \code{data} list, then these are all merged into a single table. Otherwise, it contains just the gene names and which list they were present in.
}
\section{Methods}{
\describe{

\item{\code{signature(ccGeneList = "ccGeneList", ccOptions = "ccOptions")}}{
%%  describe this method here
}
}}
\seealso{
	\code{\link{ccGeneList}} \code{\link{ccOptions}} \code{\linkS4class{mergedData}}
}
\examples{
data(ccData)
g10 <- (unique(table10$Entrez[1:100]))
g48 <- (unique(table48$Entrez[1:100]))

list10 <- list(genes=g10, universe=gUniverse, annotation="org.Hs.eg.db", data=table10[1:100,])
list48 <- list(genes=g48, universe=gUniverse, annotation="org.Hs.eg.db", data=table48[1:100,])

geneLists <- list(T10=list10, T48=list48)
geneLists <- new("ccGeneList", geneLists, ccType=c("BP","KEGG"))
ccOpts <- new("ccOptions", listNames = names(geneLists))
mergedDat <- mergeLists(geneLists,ccOpts)

list10 <- list(genes=g10, universe=gUniverse, annotation="org.Hs.eg.db")
list48 <- list(genes=g48, universe=gUniverse, annotation="org.Hs.eg.db")
geneLists <- list(T10=list10, T48=list48)
geneLists <- new("ccGeneList", geneLists, ccType=c("BP","KEGG"))
ccOpts <- new("ccOptions", listNames = names(geneLists))
mergedDat <- mergeLists(geneLists,ccOpts)
}
\keyword{methods}
