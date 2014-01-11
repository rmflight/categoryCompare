# New classes that slightly extend the original ones in Category and GOstats, enabling use of our custom function 
# and the incorporation of a false discovery rate via simulations   
setClass("HyperGParamsCC",
		contains="HyperGParams",
		representation=representation(fdr="numeric",data="data.frame"),
		prototype=prototype(annotation="GO",fdr=50,data=data.frame(0)))

setClass("GOHyperGParamsCC",
		contains=c("GOHyperGParams","HyperGParamsCC"))

setClass("KEGGHyperGParamsCC",
		contains=c("KEGGHyperGParams","HyperGParamsCC"))

setClass("HyperGResultCC",
		contains="HyperGResult",
representation=representation(minCount="numeric",fdr="numeric",fdrvalues="numeric",pvalueType="character",data="data.frame"))
		
setClass("GOHyperGResultCC",
		contains="HyperGResultCC")

setClass("KEGGHyperGResultCC",
		contains="HyperGResultCC")

# this class stores all of the various options such as which lists to compare, and the coloring, and changes colors, cssClasses,
# and compareIndx automatically, making it useful for changing options quickly and re-running the comparisons.
setClass("ccOptions",
				 representation(listNames="character",
							 compareNames="character",
							 compareIndx="list",
							 compareColors="character",
               unsaturatedColor="character",
               colorType="character",
							 cssClass="character",
							 outType="character"),
				 prototype=prototype(outType="none", colorType="solid"))

# this actually contains the actual lists of genes, gene universes, annotations, etc that the various calculations will be done on
# as part of categoryCompare. It is really just an extension of "lists", but by defining it we can set up methods
# that will check that all the various pieces are there.
setClass("ccGeneList", contains="namedList",
    representation=representation(fdr="numeric", pvalueCutoff="numeric", ccType="character",
    															testDirection="character"),
    prototype=prototype(fdr=50, pvalueCutoff=0.05, ccType="BP", testDirection="over"))
    
setClass("ccEnrichResult", contains="namedList", 
				 representation(minCount="numeric",fdr="numeric", pvalueCutoff="numeric",
				 							 pvalueType="character", categoryName="character", ontology="character",
				 							 graphType="character"))

setClass("GOccEnrichResult", contains="ccEnrichResult",
				 representation=representation(ontology="character"),
				 prototype=prototype(categoryName="GO",graphType="overlap"))

setClass("KEGGccEnrichResult", contains="ccEnrichResult", prototype=prototype(categoryName="KEGG",graphType="overlap"))
    
setClass("ccEnrichCollection", contains="namedList")

setClass("ccCompareResult", representation=representation(mainGraph="graph", subGraph="list", mainTable="data.frame", allAnnotation="list", categoryName="character", ontology="character"))

setClass("ccCompareCollection", contains="namedList")

setClass("ccSigList", representation=representation(sigID="character",categoryName="character",ontology="character",annotation="character"))

setClass("GENccEnrichResult", contains="namedList", representation=representation(categoryName="character", ontology="character", geneAnnMapping="namedList",graphType="character"), prototype=prototype(graphType="overlap") )

setClass("mergedData", contains="data.frame", representation=representation(useIDName="character"))
