# generic functions for categoryCompare
setGeneric("hyperGTestCC",
  		function(p) standardGeneric("hyperGTestCC"),
			valueClass="HyperGResultCC")

# This is the method for accessing and setting the number of replicate runs to perform
setGeneric("fdr", function(object) standardGeneric("fdr"))

setGeneric("fdr<-", function(object, value) standardGeneric("fdr<-"))

setGeneric("fdrvalues", function(object) standardGeneric("fdrvalues"))

setGeneric("fdrvalues", function(object) standardGeneric("fdrvalues"))

setGeneric("pCC", signature="object", function(object, pvalueType) standardGeneric("pCC"))

# get and modify the type of p-values to use in filtering the results
setGeneric("pvalueType", function(object) standardGeneric("pvalueType"))

setGeneric("pvalueType<-", function(object, value) standardGeneric("pvalueType<-"))

setGeneric("listNames", signature="object", function(object) standardGeneric("listNames"))

setGeneric("listNames<-", signature="object", function(object, value) standardGeneric("listNames<-"))

setGeneric("compareNames", signature="object", function(object) standardGeneric("compareNames"))

setGeneric("compareNames<-", function(object, value) standardGeneric("compareNames<-"))

setGeneric("compareIndx", function(object) standardGeneric("compareIndx"))

setGeneric("compareIndx<-", function(object, value) standardGeneric("compareIndx<-"))

setGeneric("compareColors", function(object) standardGeneric("compareColors"))

setGeneric("compareColors<-", function(object, value) standardGeneric("compareColors<-"))

setGeneric("cssClass", function(object) standardGeneric("cssClass"))

setGeneric("makeValidccOptions", function(object) standardGeneric("makeValidccOptions"))

setGeneric("makeValidccLists", function(object) standardGeneric("makeValidccLists"))

setGeneric("ccType", function(object) standardGeneric("ccType"))

setGeneric("ccType<-", function(object, value) standardGeneric("ccType<-"))

setGeneric("outType", function(object) standardGeneric("outType"))

setGeneric("outType<-", function(object, value) standardGeneric("outType<-"))

setGeneric("minCount", function(object) standardGeneric("minCount"))

setGeneric("minCount<-", function(object, value) standardGeneric("minCount<-"))

setGeneric("ccEnrich", function(ccGeneList) standardGeneric("ccEnrich"))

setGeneric("ccCompare", function(ccEnrichResult, ccOptions) standardGeneric("ccCompare"))

# setGeneric("ccCompareGO", function(GOccEnrichResult, ccOptions) standardGeneric("ccCompareGO"))

# setGeneric("ccCompareKEGG", function(KEGGccEnrichResult, ccOptions) standardGeneric("ccCompareKEGG"))

setGeneric("mainTable", function(object) standardGeneric("mainTable"))

setGeneric("allAnnotation", function(object) standardGeneric("allAnnotation"))

setGeneric("mainGraph", function(object) standardGeneric("mainGraph"))

setGeneric("ccOutCyt", function(ccCompRes,ccOpts, ...) standardGeneric("ccOutCyt"))

setGeneric("ccCompareGeneric", function(gccResult,ccOptions) standardGeneric("ccCompareGeneric"))

# setGeneric("category", function(object) standardGeneric("category"))

setGeneric("categoryName<-", function(object, value) standardGeneric("categoryName<-"))

setGeneric("sigID", function(object) standardGeneric("sigID"))

setGeneric("geneAnnMapping", function(object) standardGeneric("geneAnnMapping"))

setGeneric("graphType", function(object) standardGeneric("graphType"))

setGeneric("graphType<-", function(object, value) standardGeneric("graphType<-"))

setGeneric("mergeLists", function(ccGeneList,ccOptions,isGene=TRUE) standardGeneric("mergeLists"))

setGeneric("breakEdges", function(cwObject,cutoff,edgeAtt='weight',valDir='under',layout='force-directed') standardGeneric("breakEdges"))

setGeneric("cwReload", function(oldCW,ccOpts,...) standardGeneric("cwReload"))

setGeneric("resetColors", function(cwObj,ccOpts,...) standardGeneric("resetColors"))

setGeneric("minNodes", function(cwObj,cutoff) standardGeneric("minNodes"))

setGeneric("cytOutNodes", function(descStr, cwObj, saveObj=vector('list', 0), outImages=NULL) standardGeneric("cytOutNodes"))

setGeneric("cytOutData", function(saveObj, compareResult=NULL, mergedData=NULL, orgType="header", fileName=NULL, displayFile=FALSE) standardGeneric("cytOutData"))
