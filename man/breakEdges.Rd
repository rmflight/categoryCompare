\name{breakEdges}
\alias{breakEdges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Break Cytoscape (or graphNEL) Network Edges
}
\description{
Removes those edges in a graph network with the edge attribute that is "under" or "over" the cutoff value supplied.
}
\usage{
breakEdges(cwObject,cutoff,edgeAtt='weight',valDir='under',layout='force-directed')

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cwObject}{
The object returned by \code{\link{ccOutCyt}}
}

  \item{cutoff}{
What is the cutoff to select edges
}
  \item{edgeAtt}{
Which edge attribute should be used for deciding which edges to select
}
  \item{valDir}{
Select the edge attribute with values "under" (default) or "over" the input value
}
  \item{layout}{
What type of layout should be used after the edges are broken
}

}

\details{
When viewing annotations returned by \code{\link{ccCompare}}, unless a hierarchical view is used, the annotations are linked by edges weighted by how many genes are shared between each annotation. The default cutoff is 0.1, any edges above that are kept. This can result in a big hairy mess when viewing the annotation network in Cytoscape. To help clean up the display and encourage finding functionally related sets of annotations, it is very useful to remove edges with low values of overlap, and re-layout the network. \code{breakEdges} does just that.

The other option is to break the highly overlapping edges prior to sending the network to Cytoscape, thereby speeding up the process, as the fewer edges one has to "push" to Cytoscape the faster it will go.

}
\value{
The first returns nothing, rather it modifies the network in Cytoscape itself.

The second returns a new \code{ccCompareResult} object for use by \code{\link{ccOutCyt}}
}

\author{
Robert M Flight
}

\seealso{
\code{\link{ccOutCyt}} \code{\link{breakEdges-methods}}
}
\examples{
\dontrun{
	require(RCy3)
	g <- makeSimpleGraph()
	cw <- CytoscapeWindow('breakEdges',g)
	displayGraph(cw)
	layout(cw)
	redraw(cw)
	breakEdges(cw,'score',0,layout=NULL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
