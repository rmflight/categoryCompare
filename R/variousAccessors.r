
# various accessors for ccGeneList
setMethod("fdr", "ccGeneList", function(object) object@fdr)
setReplaceMethod("fdr", "ccGeneList", function(object, value) {
  object@fdr <- value
  object
})

# generic defined in "Category"
setMethod("pvalueCutoff", "ccGeneList", function(r) r@pvalueCutoff)

setReplaceMethod("pvalueCutoff", "ccGeneList", function(r, value) {
  r@pvalueCutoff <- value
  r
})

setMethod("ccType", "ccGeneList", function(object) object@ccType)
setReplaceMethod("ccType", "ccGeneList", function(object, value) {
  object@ccType <- value
  object
})

setMethod("listNames", "ccGeneList", function(object) (names(object)))

# generic defined in "Category"
setMethod("testDirection", "ccGeneList", function(r) r@testDirection)
setReplaceMethod("testDirection", "ccGeneList", function(r, value){
  r@testDirection <- value
  r
})

# and now some accessors for some of the results

# add in the FDR runs for HyperGParams and HyperGResult
setMethod("fdr", "HyperGParamsCC", function(object) object@fdr)
setReplaceMethod("fdr", "HyperGParamsCC", function(object, value) {
  object@fdr <- value
  object
})

setMethod("fdr", signature("HyperGResultCC"), function(object) object@fdr)

setMethod("fdrvalues", "HyperGResultCC", function(object) object@fdrvalues)


# give the ability to modify the pvalueCutoff in HyperGResult objects (makes it easier to see changes using "show")
setReplaceMethod("pvalueCutoff", "HyperGResultCC", function(r, value) {
  r@pvalueCutoff <- value
  r
})

# want a method for pvalueType for the HyperGResultCC objects
setMethod("pvalueType", "HyperGResultCC", function(object) object@pvalueType)
setReplaceMethod("pvalueType", "HyperGResultCC", function(object, value) {
  allowedType <- c("fdr","pval")
  if (value %in% allowedType){
    object@pvalueType <- value
  } else {
    stop("pvalueType must be one of either: ", paste(allowedType, collapse="; "), call.=FALSE)
  }
  object
})

# these are used to actually get our results our of the HyperGResult objects themselves
isResult <- function(result, pvalue, pvalueType, minCount=NULL){
  pvals <- pCC(result, pvalueType)
	wanted <- is.finite(pvals) & pvals <= pvalue
	if (!is.null(minCount)){
		gCounts <- geneCounts(result)
		hasMinSize <- gCounts >= minCount
		wanted <- wanted & hasMinSize
	}
	wanted
}

# we should have a method that returns the annotated genes from our object as well, so we don't have to do it over and over again in the main software package 
# ** This actually is already in the geneIdsByCategory function in Category

# this is the pvalues that can be selected using either by pvalueType of "fdr" or "pval".
setMethod("pCC", signature(object="HyperGResultCC"),
		function(object, pvalueType) {
			if ((missing(pvalueType)) || (is.null(pvalueType))){
				pvalueType <- pvalueType(object)
			}
			if (pvalueType == 'fdr'){
				fdrvalues(object)
			} else { pvalues(object) }
		})

# the replace methods essentially propogate down to the objects contained within
setMethod("fdr", "ccEnrichResult", function(object) object@fdr)

setMethod("pvalueCutoff", "ccEnrichResult", function(r) r@pvalueCutoff)
setReplaceMethod("pvalueCutoff", "ccEnrichResult", function(r, value) {
  r@pvalueCutoff <- value
  
  nSub <- length(r)
  for (iSub in 1:nSub){
    pvalueCutoff(r[[iSub]]) <- value
  }
  r
})

setReplaceMethod("pvalueCutoff", "ccEnrichCollection", function(r, value) {
  nSub <- length(r)
  for (iSub in 1:nSub){
    pvalueCutoff(r[[iSub]]) <- value
  }
  r
})

setReplaceMethod("pvalueType", "ccEnrichResult", function(object, value) {
  nSub <- length(object)
  for (iSub in 1:nSub){
    pvalueType(object[[iSub]]) <- value
  }
  object@pvalueType <- value
  object
})

setReplaceMethod("pvalueType", "ccEnrichCollection", function(object, value) {
  nSub <- length(object)
  for (iSub in 1:nSub){
    pvalueType(object[[iSub]]) <- value
  }
  object
})

setMethod("pvalueType","ccEnrichResult", function(object) object@pvalueType)

## should set a method for minCount as a piece of data that we might want to use
setMethod("minCount", "ccEnrichResult", function(object) object@minCount)

setMethod("minCount", "HyperGResultCC", function(object) object@minCount)

setReplaceMethod("minCount", "ccEnrichResult", function(object, value) {
  nSub <- length(object)
  for (iSub in 1:nSub){
    minCount(object[[iSub]]) <- value
  }
  object@minCount <- value
  object
})

setReplaceMethod("minCount", "ccEnrichCollection", function(object, value) {
  nSub <- length(object)
  for (iSub in 1:nSub){
    minCount(object[[iSub]]) <- value
  }
  object
})

setReplaceMethod("minCount", "HyperGResultCC", function(object, value) {
	object@minCount <- value
	object
	})


setMethod("graphType", "ccEnrichResult", function(object) {
	object@graphType
})

setReplaceMethod("graphType", "ccEnrichResult", function(object, value) {
	object@graphType <- value
  object
})

setReplaceMethod("graphType", "ccEnrichCollection", function(object, value) {
	nSub <- length(object)
  for (iSub in 1:nSub){
  	graphType(object[[iSub]]) <- value
  }
  object
})

## categoryName for some different objects
setMethod("categoryName", "ccEnrichResult", function(r){
	r@categoryName
})

setReplaceMethod("categoryName", "ccEnrichResult", function(object, value){
	object@categoryName <- value
	object
})

setMethod("categoryName", "ccCompareResult", function(r){
	r@categoryName
})


### ccSigList objects
setMethod("sigID", "ccSigList", function(object) {
	object@sigID
})

setMethod("categoryName", "ccSigList", function(r) {
	r@categoryName
})

setReplaceMethod("categoryName", "ccSigList", function(object, value) {
	object@categoryName <- value
	object
})

setMethod("ontology", "ccSigList", function(object) {
	object@ontology
})

setMethod("annotation", "ccSigList", function(object) {
	object@annotation
})

### GENccEnrichResult objects
setMethod("categoryName", "GENccEnrichResult", function(r) {
	r@categoryName
})

setReplaceMethod("categoryName", "GENccEnrichResult", function(object, value) {
	object@categoryName <- value
	object
})

setMethod("ontology", "GENccEnrichResult", function(object) {
	object@ontology
})

setMethod("geneAnnMapping", "GENccEnrichResult", function(object) {
	object@geneAnnMapping
})

setMethod("graphType", "GENccEnrichResult", function(object) {
	object@graphType
})

setReplaceMethod("graphType", "GENccEnrichResult", function(object, value) {
	object@graphType <- value
  object
})

# Needed for accessing sub-pieces of the ccEnrichResult objects to allow a person to do only sub-comparisons
# because it is really just an extension of lists, we lose the important slots if we don't specify a method for it.
setMethod("[", "ccEnrichResult", function(x, i){
	currClass <- class(x)	
	if (class(i) == "character"){
		i <- match(i,names(x),nomatch=0)
	}
	xRet <- new("namedList")
	xRet@.Data <- x@.Data[i]
	names(xRet) <- names(x)[i]
	xRet <- new(currClass,xRet,
							minCount=x@minCount,
							fdr=x@fdr,
							pvalueCutoff=x@pvalueCutoff,
							pvalueType=x@pvalueType,
							categoryName=x@categoryName,
							ontology=x@ontology,
							graphType=x@graphType)
	return(xRet)
})

setMethod("[", "GENccEnrichResult", function(x, i){
	currClass <- class(x)	
	if (class(i) == "character"){
		i <- match(i,names(x),nomatch=0)
	}
	xRet <- new("namedList")
	xRet@.Data <- x@.Data[i]
	names(xRet) <- names(x)[i]
	xRet <- new(currClass,xRet,
							categoryName=x@categoryName,
							ontology=x@ontology,
							geneAnnMapping=x@geneAnnMapping,
							graphType=x@graphType)
	return(xRet)
})