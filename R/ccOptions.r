# methods for ccOptions
# this helps us to call a function to make sure that our object is created properly.
setMethod("initialize", "ccOptions",
  function(.Object, ...){
  .Object <- callNextMethod()
  .makeValidccOptions(.Object)
})

.makeValidccOptions <- function(object){
 # check and set up the ccOptions object so that it is useable by later functions
 # biggest thing is figuring stuff out if the user has put "all" for compareNames or supplied actual stuff 
 # for compareNames and compareColor. Also need to figure out cssClass (not set by user) so we can color our HTML tables
 # correctly.
 currNames <- listNames(object)
 if (length(currNames) == 0){
   stop("listNames must be supplied!", call.=FALSE)
 }
 if (class(currNames) != "character"){
   stop("listNames must be a vector of strings to denote which lists will be used in the comparisons.", call.=FALSE)
 }
 
 currCompare <- compareNames(object)
 if (length(currNames) == 0){
   warning('No compareNames specified, assuming "all"')
   compareNames(object) <- 'all'
 } else {
   compareNames(object) <- currCompare
 }
 currColors <- compareColors(object)
 compareColors(object) <- currColors
 outType(object) <- outType(object) # this should just set defaults
 return(object)
  
}

# setMethod("makeValidccOptions", "ccOptions", .makeValidccOptions)

setMethod("listNames", "ccOptions", function(object) object@listNames)
setReplaceMethod("listNames", "ccOptions", function(object, value) {
  object@listNames <- value
	object
})

setMethod("compareNames", "ccOptions", function(object) object@compareNames)
setReplaceMethod("compareNames", "ccOptions", function(object, value) .replaceCompareNames(object, value))

.replaceCompareNames <- function(object, value){
  currNames <- listNames(object)
  if (colorType(object) == "solid"){
    if ((length(value)==0) || (tolower(value[1]) == "all")){
      compData <- .compData(currNames)
    } else {
      compData <- .compData(currNames,value)
    }
    # can't use the functions because otherwise we will just enter an infinite recursion, which is not good
    object@compareNames <- compData$name
    object@compareIndx <- compData$indx
    
    object@cssClass <- .classGen(compData$name)
  } else {
    object@compareNames <- currNames
    tmpIndx <- seq(1, length(currNames))
    names(tmpIndx) <- currNames
    object@compareIndx <- as.list(tmpIndx)
  }
    
  object
}

setMethod("compareIndx", "ccOptions", function(object) object@compareIndx)
setReplaceMethod("compareIndx", "ccOptions", function(object, value) {
  object@compareIndx <- value
  object
})

# actually does the checking of compareNames against the listNames, and generates new comparisons and indices if required
.compData <- function(tmpNames,compName=NULL){
  
	# check if all of our listNames are in the compNames
  
	# either we have to figure out all of the comparisons, or search through to find the indices
	if (is.null(compName)){
		compIndx <- vector("list", 0)
		compName <- vector("character", 0)
		nList <- length(tmpNames)
		listSeq <- seq(1,nList)
		
		for (iList in 1:nList){
			compIndx <- append(compIndx,combn(listSeq, iList, simplify=F))
		}
		compName <- sapply(compIndx, function(x) {
			paste(tmpNames[x], collapse=',')
		})
	} else {
		nComp <- length(compName)
		splitName <- strsplit(compName, ',')
		compIndx <- lapply(splitName, function(x) {
			match(x,tmpNames)})
	}
	names(compIndx) <- compName
	compData <- list(indx=compIndx, name=compName)
}

setMethod("cssClass", "ccOptions", function(object) object@cssClass)

# turn out classnames to use in the CSS to define colors based on which grouping genes belong to
# basically just strip out the commas to create long string variable names
.classGen <- function(orgNames)
{
  nName <- length(orgNames)
	newNames <- sapply(orgNames, function(x){
		tmpName <- strsplit(x,"[[:punct:]]")
		paste(tmpName[[1]],collapse='')
	})
	names(newNames) <- orgNames
	newNames
}

setMethod("compareColors", "ccOptions", function(object) object@compareColors)
setReplaceMethod("compareColors", "ccOptions", function(object,value) {
  .colorData(object,value)
  })

# takes potential list of colors, these can be hex, rgb, or colornames
# checks that we have enouch colors
.colorData <- function(object,compColor=NULL){
  
  # the comparisons we want to make
  compNames <- compareNames(object)
  
  # what kind of colors did we pass in? Matrix of RGB, or character vector
  colorClass <- class(compColor)
  
  nComp <- length(compNames)
	allColors <- colors() # these are the colors (text strings) that R knows about
	makeColor <- TRUE
		
  # if nothing, then figure out the colors  
	if (is.null(compColor)){
		makeColor <- TRUE
	} else if (length(grep('#',compColor)) >= nComp){  
	# check for hex colors ( and assume they are properly formatted)
		makeColor <- FALSE
    compColor <- compColor[1:nComp]
	}
	# if the names in compColor
	else if (sum(tolower(compColor) %in% allColors) >= nComp){
		makeColor <- FALSE
		compColor <- compColor[1:nComp]
	}

	if (makeColor){	
		endCol <- 315 * (nComp - 1) / nComp
		compColor <- rainbow_hcl(nComp, c=100, start=0, end=endCol)
	}
	
	names(compColor) <- compNames
	object@compareColors <- compColor
  
  if (colorType(object) == "pie"){
    unSatColor <- desaturate(compColor)
    names(unSatColor) <- compNames
    object@unsaturatedColor <- unSatColor
  }
  object

}

setMethod("outType", "ccOptions", function(object) object@outType)
setReplaceMethod("outType", "ccOptions", function(object,value) {
  .outType(object,value)
})

.outType <- function(object,value){
  validTypes <- c('html','text','rcytoscape')
  validNone <- 'none'
  validAll <- c(validTypes,validNone)
  nVal <- length(validAll)
  if (missing(value) | is.null(value)){
    value <- 'none'
    message("No valid out types supplied, setting to none")
  } else {
    value <- tolower(value)
    keepVal <- value %in% validAll
    isNone <- value %in% validNone
    if (sum(keepVal) == 0){
      value <- 'none'
      message("No valid out types supplied, setting to none")
    } else {
    	if ((sum(keepVal) > 1) & (sum(isNone) > 0)){
    		keepVal[isNone] <- FALSE
    	}
      value <- value[keepVal]
    }
  }
  object@outType <- value
  object
}

setMethod("colorType", "ccOptions", function(object) object@colorType)