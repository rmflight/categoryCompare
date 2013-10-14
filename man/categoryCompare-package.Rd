\name{categoryCompare-package}
\alias{categoryCompare-package}
\alias{categoryCompare}
\docType{package}
\title{
Meta-analysis of high-throughput experiments using feature annotations
}
\description{
Calculates significant annotations (categories) in each of two (or more) feature (i.e. gene) lists, determines the overlap between the annotations, and returns graphical and tabular data about the significant annotations and which combinations of feature lists the annotations were found to be significant. Interactive exploration is facilitated through the use of RCytoscape (heavily suggested).
}
\details{
\tabular{ll}{
Package: \tab categoryCompare\cr
Version: \tab 0.99.1\cr
License: \tab GPL-2\cr
Depends: \tab Biobase (>= 1.15.29), AnnotationDbi (>= 0.1.15), Category\cr
Suggests: \tab methods, GSEABase, hwriter, colorspace, graph, GO.db, KEGG.db, estrogen, org.Hs.eg.db, hgu95av2.db\cr
Imports: \tab Biobase (>= 1.15.29), AnnotationDbi (>= 0.1.15), hwriter, GSEABase, Category (>= 2.21.2), GOstats, annotate, colorspace, graph, RCytoscape (>= 1.5.11)\cr
LazyLoad: \tab yes\cr
biocViews: \tab Bioinformatics, Annotation, GO, MultipleComparisons, Pathways, GeneExpression\cr
SystemRequirements: \tab Cytoscape (>= 2.8.0) (if used for visualization of results, heavily suggested), CytoscapeRPC plugin (>= 1.8)\cr
TODO: \tab Text and HTML output without graphs.\cr
Built: \tab R 2.15.0; ; 2012-03-15 18:42:40 UTC; windows\cr
}

Index:
\preformatted{
GENccEnrichResult-class
                        Class '"GENccEnrichResult"'
HyperGParamsCC-class    Class "HyperGParamsCC"
HyperGResultCC-class    Class "HyperGResultCC"
breakEdges              Break Cytoscape (or graphNEL) Network Edges
breakEdges-methods      Methods for Function 'breakEdges' in Package
                        'categoryCompare'
categoryCompare-package
                        Meta-analysis of high-throughput experiments
                        using feature annotations
ccCompare-methods       Comparison of enriched annotations
ccCompareCollection-class
                        Class '"ccCompareCollection"'
ccCompareGeneric-methods
                        Methods for Function 'ccCompareGeneric' in
                        Package 'categoryCompare'
ccCompareResult-class   Class '"ccCompareResult"'
ccData                  Test data for 'categoryCompare'
ccEnrich-method         Perform annotation enrichment for multiple gene
                        lists
ccEnrichCollection-class
                        Class "ccEnrichCollection"
ccEnrichResult-class    Class "ccEnrichResult"
ccGeneList-class        Class "ccGeneList"
ccOptions-class         Class "ccOptions"
ccOutCyt-methods        Methods for Function 'ccOutCyt' in Package
                        'categoryCompare'
ccSigList-class         Class '"ccSigList"'
cwReload-methods        Methods for Function 'cwReload' in Package
                        'categoryCompare'
cytOutData-methods      Methods for Function 'cytOutData'
cytOutNodes-methods     Methods for Function 'cytOutNodes'
fdr                     Number of FDR runs to perform
getGeneSymbol           Entrez to name, symbol, GO and path conversion,
                        as well as general ID to ID conversion.
graphType-methods       graphType
hyperGTestCC            Hypergeometric testing with false discovery
                        rate
listNames               listNames
mergeLists-methods      Function 'mergeLists' in Package
                        'categoryCompare'
mergedData-class        Class '"mergedData"'
minCount                minCount
minNodes                Delete nodes with less than a certain number of
                        genes annotated
pvalueType              Type of p-values to return from object
resetColors-methods     resetColors
show-methods            Methods for Function show in Package
                        'categoryCompare'
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{categoryCompare_vignette} \tab categoryCompare: High-throughput data meta-analysis using gene annotations (source)\cr
}
}
\author{
Robert M. Flight <rflight79@gmail.com>

Maintainer: Robert M. Flight <rflight79@gmail.com>
}

\keyword{ package }