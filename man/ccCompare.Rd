\name{ccCompare-methods}
\alias{ccCompare}
\alias{ccCompare-methods}
\alias{ccCompare,ccEnrichCollection,ccOptions-method}
\alias{ccCompare,GOccEnrichResult,ccOptions-method}
\alias{ccCompare,KEGGccEnrichResult,ccOptions-method}
\alias{ccCompare,GENccEnrichResult,ccOptions-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Comparison of enriched annotations
}
\description{
Takes the results from \code{\link{ccEnrich}} and compares the enriched annotations based on the settings previously set in \code{\link{ccOptions}}. Returns a \code{ccCompareResult} or \code{ccCompareCollection} object, see Details.
}
\usage{
ccCompare(ccEnrichResult, ccOptions)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ccEnrichResult}{
The enriched annotations collection returned from \code{ccEnrich}. This can be the \code{ccEnrichCollection}, \code{GOccEnrichResult}, or \code{KEGGccEnrichResult}
}
  \item{ccOptions}{
A \code{\link{ccOptions}} object that will determine which lists to actually compare against each other. See details below.
}

}
\details{
Based on the enrichments found for each gene list, we now want to compare the annotations between lists. \code{ccCompare} accesses the annotations for each enrichment performed for each list, and makes the comparisons defined in \code{ccOptions}.

}
\value{
\code{ccCompare} generates both a graph of the comparisons (to show how the categories are linked to each list and each other) and tabular output. The tabular output is a data frame, with \code{ID} for each term that was considered as a candidate annotation for each list, as well as a long description (\code{Desc}) of what the term is, and then membership and statistics from each gene list.

For each type of comparison (\code{GO}, \code{KEGG}, etc) a \code{ccCompareResult} is generated, with the following slots:
\item{mainGraph }{Annotations arranged as a graph}
\item{mainTable }{The tabular results from all enrichment calculations combined into one}
\item{allAnnotation }{A list of lists, where each entry is the annotation identifier, then a list for each comparison, with the genes that are annotated to that term that also belong to each list}

The default is to generate an overlap graph for GO and KEGG, where the overlap is a measure of the similarity of the features (genes) annotated to each annotation term (based on a formula from \code{EnrichmentMap}). Optionally for GO, one can generate a hierarchical layout where the parent GO terms of the significant terms will also be included in the graph, with term origin saved in the node annotation (see example below to do this).

 Only those terms with more than 10 and less than 500 annotated genes (according to the GO annotation file) are included.

When using weighted overlap graphs and \pkg{RCytoscape} for viewing, it is recommended to use \code{breakEdges} and \code{minNodes} to remove edges with low weights and nodes with only a few genes from the dataset annotated to them.

}

\author{
Robert M Flight
}

\seealso{
\code{\linkS4class{ccCompareResult}} \code{\linkS4class{ccCompareCollection}} \code{\link{ccOutCyt}} \code{\link{breakEdges}} \code{\link{outType}} \code{\link{ccEnrich}}}

\examples{
\dontrun{
require(GO.db)
require(KEGG.db)
require(org.Hs.eg.db)
}
data(ccData)

# note that enrichLists is generated from ccEnrich
# ccResults <- ccCompare(enrichLists,ccOpts)
ccResults

# use the GO hierarchy tree
graphType(enrichLists$CC) <- "hierarchical"
# ccResultsHier <- ccCompare(enrichLists$CC,ccOpts)
ccResultsHier


}
