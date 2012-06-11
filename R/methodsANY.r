# these are replacement functions for many of the Category methods to work on our ANYHyperGParams objects
.makeValidParams <- function(object){
	sel <- geneIds(object)
	if (is.list(sel)) {
		warning("converting geneIds from list to atomic vector via unlist")
		sel <- unlist(sel)
	}
	if (any(duplicated(sel))) {
		warning("removing duplicate IDs in geneIds")
		sel <- unique(sel)
	}
	geneIds(object) <- sel
	
	univ <- universeGeneIds(object)
	if (length(univ)) {
		if (is.list(univ)){
			warning("converting univ from list to atomic vector via unlist")
			univ <- unlist(univ)
		}
		if (typeof(sel) != typeof(univ)){
			stop(paste("geneIds and universeGeneIds must have the same mode\n",
								 "geneIds:", typeof(sel), "\n",
								 "universeGeneIds:", typeof(univ)), .Call=FALSE)
		}
		if (any(duplicated(univ))){
			warning("removing duplicate IDs in universeGeneIds")
			univ <- unique(univ)
		}
		universeGeneIds(object) <- univ
		if (!all(sel %in% univ)){
			warning("removing geneIds not in universeGeneIds")
			sel <- intersect(sel, univ)
			if (!length(sel)){
				stop("no geneIds in universeGeneIds", .Call=FALSE)
			}
			geneIds(object) <- sel
		}
	}
	pv <- pvalueCutoff(object)
	if (pv > 1 || pv < 0){
		stop("invalid pvalueCutoff, must be between 0 and 1", .Call=FALSE)
	}
	return(object)
}

setMethod("makeValidParams", "ANYHyperGParamsCC", .makeValidParams)

setMethod("geneIds", "ANYHyperGParamsCC", function(object, ...) object@geneIds)
setReplaceMethod("geneIds", "ANYHyperGParamsCC", function(object, value) {
    object@geneIds <- value
    object
})


setMethod("testDirection", "ANYHyperGParamsCC", function(r) r@testDirection)
setReplaceMethod("testDirection", "ANYHyperGParamsCC", function(r, value) {
    r@testDirection <- value
    r
})


setMethod("categoryName", "ANYHyperGParamsCC", function(object) object@categoryName)
setReplaceMethod("categoryName", "ANYHyperGParamsCC", function(object, value) {
    object@categoryName <- value
    object
})

setMethod("annotation", "ANYHyperGParamsCC", function(object) object@annotation)
setReplaceMethod("annotation", c("ANYHyperGParamsCC", "list"),
                 function(object, value) {
                   object@annotation <- value
                   object
                 })

setMethod("pvalueCutoff", "ANYHyperGParamsCC", function(r) r@pvalueCutoff)
setReplaceMethod("pvalueCutoff", "ANYHyperGParamsCC", function(r, value) {
    r@pvalueCutoff <- value
    r
})

setMethod("universeGeneIds", "ANYHyperGParamsCC", function(r) r@universeGeneIds)
setMethod("universeGeneIds<-", "ANYHyperGParamsCC", function(r, value) {
    r@universeGeneIds <- value
    r
})

setMethod("fdr", "ANYHyperGParamsCC", function(object) object@fdr)
setReplaceMethod("fdr", "ANYHyperGParamsCC", function(object, value) {
	object@fdr <- value
	object
})

setMethod("conditional", "ANYHyperGParamsCC", function(r) r@conditional)

setMethod("organism", "ANYHyperGParamsCC", function(object) object@organism)
### now methods for setting up the universe and categoryToEntrez stuff
setMethod("universeBuilder", signature(p="ANYHyperGParamsCC"),
          function(p) {
                getUniverseANY(p)
          })
# assumes you are passing a list object with categories as the list with each containing a vector of identifiers annotated to that category    
getUniverseANY <- function(object){
	tmpCategory <- annotation(object)
	univ <- universeGeneIds(object)
	id2Cat <- reverseSplit(tmpCategory)
	id2CatLen <- sapply(id2Cat, length)
	id2Cat <- id2Cat[id2CatLen > 0]
	univ <- intersect(univ,names(id2Cat))
	if (!length(univ)){
		stop("no universeGeneIds annotated to any Category", .Call=FALSE)
		
	} else {
		return(univ)
	}
}
    
setMethod("categoryToEntrezBuilder",
          signature(p="ANYHyperGParamsCC"),
          function(p) {
              getANYidMap(p)
          })
# returns only those categories that have something in the Universe mapped to it, and removes those id's in the map that are not in the universe
getANYidMap <- function(object){
	univ <- universeGeneIds(object)
	id2Cat <- reverseSplit(annotation(object))
	idKeep <- intersect(univ, names(id2Cat))
	id2Cat <- id2Cat[idKeep]
	cat2ID <- reverseSplit(id2Cat)
	if (!length(cat2ID)){
		stop("no categories annotated to genes in the universe")
	} else{ 
		return(cat2ID)
	}
}