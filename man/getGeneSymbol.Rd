\name{getGeneSymbol}
\alias{getGeneSymbol}
\alias{getGeneName}
\alias{getGO2ALLEGS}
\alias{getPATH2EG}
\alias{getAnnotation}
\title{
Entrez to name, symbol, GO and path conversion, as well as general ID to ID conversion.
}
\description{
Get different attributes for the Entrez gene Ids
}
\usage{
getGeneSymbol(id, annPackage)
getGeneName(id, annPackage)
getGO2ALLEGS(id, annPackage)
getPATH2EG(id, annPackage)
getAnnotation(id, annPackage, mapID, doUnlist=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{The IDs one wants to get information for.
}
  \item{annPackage}{
Which annotation package to use.
}
	\item{mapID}{Which mapping to use}
	\item{doUnlist}{should the results be unlisted or not?}
}
\details{
  The type of ID will change depending on the function. For \code{getGene...} the ID should be Entrez IDs. For \code{getGO2ALLEGS} Gene Ontology IDs should be used, and for \code{getPATH2EG} KEGG pathways IDs should be used. For \code{getAnnotation}, any ID can be used.
}
\value{
Returns the requested information.
}
\author{
Robert M Flight
}
\note{
These functions are generally called internally for mapping between genes and various objects.
}
