# Some useful utilities. Note, much of the ideas for this are borrowed from the oligoClasses package in Bioconductor

# creates a bar across the screen, useful for sending messages to the user.
getBar <- function(width=getOption("width"))
	paste(rep("=", width), collapse="")

# checks if a package is loaded
isPackageLoaded <- function(pkg){
	stopifnot(is.character(pkg))
	pkg <- paste("package:", pkg, sep="")
	pkg %in% search()
}

visStatus <- function(verbose=FALSE){
	isCy <- isPackageLoaded("RCy3")
	if (verbose){
		message(getBar())
		message("Interactive visualization support for categoryCompare: ")
		message("")
		message("Using Cytoscape:")
		if (isCy){
			message("Fully Enabled: access to all ")
		} else {
			message("Disabled")
			message("  - Load 'RCy3'")
		}
		message("")
		message(getBar())
	}
	isVis <- list(cyt=isCy)
	return(isVis)
}

annStatus <- function(verbose=FALSE){
	isGODB <- isPackageLoaded("GO.db")
	isKEGGDB <- isPackageLoaded("KEGG.db")
# 	isGOSTATS <- isPackageLoaded("GOstats")

	if (verbose){
		message(getBar())
		message("Possible annotations using currently loaded packages: ")

		message("Gene Ontology: ", appendLF=FALSE)
		if (isGODB){
			message("Enabled")
		} else {
			message("Disabled")
			message("	- Load 'GO.db'")
		}
		message("")
		message("KEGG: ", appendLF=FALSE)
		if (isKEGGDB){
			message("Enabled")
		} else {
			message("Disabled")
			if (!isKEGGDB){
				message("	- Load 'KEGG.db'")
			}
		}
		message("")
	}
	isAnn <- list(godb=isGODB,keggdb=isKEGGDB)
	return(isAnn)
}

getAnnotation <- function(id,annPackage,mapID,doUnlist=TRUE){
	if (!is.null(mapID)){
		annEnv <- getAnnMap(mapID, annPackage, load=TRUE)
	  tmpDat <- mget(id, annEnv, ifnotfound=NA)
	  if (doUnlist){
	  	tmpDat <- unlist(tmpDat, use.names=FALSE)
	  }
	} else {
		tmpDat <- vector("character", length(id))
	}
  return(tmpDat)
}

getGeneSymbol <- function(id,annPackage){
  annEnv <- getAnnMap("SYMBOL", annPackage, load=TRUE)
  unlist(mget(id, annEnv, ifnotfound=NA),use.names=FALSE)
}

getGeneName <- function(id,annPackage){
  annEnv <- getAnnMap("GENENAME", annPackage, load=TRUE)
  unlist(mget(id, annEnv, ifnotfound=NA),use.names=FALSE)
}

# this is here if we ever need it
getEntrez <- function(id,annPackage){
	annEnv <- getAnnMap("ENTREZID", annPackage, load=TRUE)
  unlist(mget(id, annEnv, ifnotfound=NA), use.names=FALSE)
}

getGO2ALLEGS <- function(id,annPackage){
  annEnv <- getAnnMap("GO2ALLEGS", annPackage, load=TRUE)
  tmpTerm <- (mget(id, annEnv, ifnotfound=NA))
  lapply(tmpTerm,unique) # get rid of repeats in each one
}

getPATH2EG <- function(id,annPackage){
  annEnv <- getAnnMap("PATH2EG", annPackage, load=TRUE)
  mget(id, annEnv, ifnotfound=NA)
}
