## ccGeneList-Methods
# define all of the methods that will act on our ccGeneList object.
# the formal use of an object is really just to make it easier to verify that everything 
# that is in our data is actually in a format that is useful, and we do this when the object is created so that the user
# does not get frustrated when the actual software does not run because the data is not in the correct format.

# this helps us to call a function to make sure that our object is created properly. Remember, it is simply a list of lists
setMethod("initialize", "ccGeneList",
  function(.Object, ...){
  .Object <- callNextMethod()
  makeValidccLists(.Object)
})



.makeValidccLists <- function(object){
  reqFields <- c(genes="character",universe="character",annotation="character")
  nReq <- length(reqFields)
  supFields <- list(data=list(class="data.frame",pos=c("ID","ENTREZ")))
  
  allowedTypes <- c("BP","MF","CC","KEGG","SPIA") # not using SPIA yet, but will 
 # check and set up the ccGeneList object so that it is useable by later functions
 # want to check that the minimum required fields are there (genes, gUniverse, annotation), and check that any other
 # pieces are also the right type and have the required bits of information. ## Note this will probably evolve over time.
 #browser(expr=TRUE)
 
 ccTypes <- object@ccType
   
 nLists <- length(object)
 currNames <- names(object) # these will always be character, it is how "names" returns them. 
 if (nLists == 0){
   stop("Empty ccGeneList. Please supply new-lists.", call.=FALSE)
 }
 
 # check that the names are not just numbers, these will not work well
 numNames <- as.numeric(currNames)
 if (sum(is.na(numNames)) == 0){
   warning('Numeric names detected, adding "V" to all object names!', call.=FALSE)
   currNames <- paste('V', currNames, sep="")
   names(object) <- currNames
 }
 
 # also check if the names begin with numbers
 numNames <- grep('^[[:digit:]]',currNames)
 if (length(numNames) > 0){
   stop('List names are not allowed to start with numbers. Please rename the lists!', call.=FALSE)
 }
 
 # run a check on each of the lists to make sure the required fields are there
 reqNames <- names(reqFields)
 reqClass <- unlist(reqFields)
 
	
 # check everything out, and stop if there is an error.
 sapply(object, function(x){
   #browser(expr=TRUE)
   subNames <- names(x)
   misNames <- reqNames[!(reqNames %in% subNames)]
   
   # check the sub-lists that are required
   if (length(misNames) > 0){
     stopStr <- paste('Missing required entry ', misNames, sep="", collapse="\n")
     stop(stopStr, call.=FALSE)
   }
   
   matchReq <- match(reqNames,subNames, nomatch=0)
   nChar <- sum(match(reqClass,(sapply(x[matchReq],"class")), nomatch=0) != 0)
   if (nReq != nChar){
     stop('Data type mismatch. Please review the documentation.', call.=FALSE)
   }
   
   # and now check for optional "data" that will be kept.
   isSup <- match(names(supFields),subNames, nomatch=0)
   isSup <- isSup[!(isSup == 0)]
   nSup <- length(isSup)
   
   if (nSup > 0){
     for (iSup in 1:nSup){
       # check the class
       if (class(x[[isSup[iSup]]]) != supFields[[iSup]]$class){
         stop('Data type mismatch.', call.=FALSE)
       }
		 	 if (supFields[[iSup]]$class == "data.frame"){
		 	 	objNames <- names(x[[isSup[iSup]]])
				posNames <- supFields[[iSup]]$pos
				if (sum(tolower(objNames) %in% tolower(posNames)) == 0){
					stop(paste('Required information missing from the extra entry ', names(supFields)[iSup], collapse=""))
				}
		 	 }
     }
   }
 })
 
 isAllowed <- ccTypes %in% allowedTypes
 nNot <- sum(!isAllowed)
 if (nNot > 0){
   notAllowed <- ccTypes[!isAllowed]
   notAllowed <- paste(notAllowed, collapse="; ")
   warnStr <- paste("The following ccTypes are not recognized, and will be excluded: ", notAllowed, sep="")
   warning(warnStr)
   object@ccType <- ccTypes[isAllowed]
 }
  
 return(object)
}
  
setMethod("makeValidccLists", "ccGeneList", .makeValidccLists)

	 
## merge data in a ccGeneList object
# isGene is used to tell us that yes, we should have the symbol, and name. If it is FALSE, then we won't bother 
# trying to future proof to be able to use things like metabolites
.mergeLists <- function(ccGeneList,ccOptions,isGene=TRUE){
	reqID <- c("entrez","symbol","name") # these are required in the final output for genes
  names(reqID) <- c("ENTREZID","SYMBOL","GENENAME")
 	optID <- "id" # if this is there, great, but doesn't have to be -> if it is there and the others aren't then we need to use it to get the others'
	geneListNames <- names(ccGeneList)
	nList <- length(ccGeneList)
  ccComps <- compareNames(ccOptions)
 	
  listAnn <- (sapply(ccGeneList,function(x) x$annotation))
  difAnn <- FALSE
  if (length(unique(listAnn)) > 1){
  	difAnn <- TRUE
  }
  
 	allDat <- data.frame(ZZ=0)
 	geneLists <- vector('list', nList)
 	names(geneLists) <- geneListNames
 	allDatUseID <- vector('character', nList) # which have the IDs from the genelists?
  for (iList in 1:nList){
  	geneLists[[iList]] <- ccGeneList[[iList]]$genes
  	# check if our primary IDs will be Entrez IDs or probe IDs
  	isEntrez <- length(grep("org",ccGeneList[[iList]]$annotation)) > 0
  	if (isEntrez && !isGene){
  		warning("Assuming Entrez IDs based on annotation, but 'isGene' supplied as FALSE. Changing to TRUE.\n
  						Change the annotation or 'isGene' to suppress this warning.")
  		isGene <- TRUE
  	}
	 	tmpID <- ccGeneList[[iList]]$genes
	  if ("data" %in% names(ccGeneList[[iList]])){
	  	# figure out what column matches our geneList
	  	tmpDat <- ccGeneList[[iList]]$data
			colClass <- sapply(tmpDat, function(x) class(x))
			potCol <- names(colClass[colClass %in% "character"])
			matchCol <- vector("integer",length(potCol))
			names(matchCol) <- potCol	
			for (iCol in 1:length(potCol)){
				matchCol[iCol] <- sum(tmpID %in% tmpDat[[potCol[iCol]]])
			}
			matchCol <- matchCol[matchCol > 0]
		  if (unique(matchCol) == 0){
		  	stop("None of the data columns match the gene list identifiers!", call.=FALSE)
		  }
			useCol <- names(matchCol)[which.max(matchCol)]
			# check if any of our identifiers are not in the data table, and add them if necessary
			misID <- tmpID[!(tmpID %in% tmpDat[,useCol])]
			if (length(misID) > 0){
				emptyDat <- tmpDat[0,]
				emptyDat[1:length(misID),] <- NA
				emptyDat[,useCol] <- tmpDat[,useCol]
				tmpDat <- rbind(tmpDat,emptyDat)
			}
	  } else {
	  	
	  	if (isEntrez){
	  		tmpDat <- data.frame(Entrez=tmpID,stringsAsFactors=FALSE)
	  	} else { 
	  		tmpDat <- data.frame(ID=tmpID,stringsAsFactors=FALSE)
		 	}
	  }
	  # check whether we have the Entrez, Symbol, and Name in our data table, but only if we know that these are genes
	  if (isGene){
	  	misReq <- reqID[!(tolower(reqID) %in% tolower(names(tmpDat)))]
			if (isEntrez){
				useID <- names(tmpDat)[tolower(names(tmpDat)) %in% "entrez"]
			} else {
				useID <- names(tmpDat)[tolower(names(tmpDat)) %in% "id"]
				tmpDat[,'entrez'] <- getAnnotation(tmpDat[,useID],listAnn[iList],"ENTREZID")
		  }
			for (iMis in 1:length(misReq)){
				tmpDat[,misReq[iMis]] <- getAnnotation(tmpDat[,useID],listAnn[iList],names(misReq)[iMis])
			}
	  } else {
	  	useID <- names(tmpDat)[tolower(names(tmpDat)) %in% "id"]
	  }
	 	allDatUseID[iList] <- useID
	 	# now merge the tables together into a single object that can be used by the user, or also used by RCytoscape to output which genes and their expression are tied to which annotation
	 	if (difAnn){
	 		names(tmpDat) <- paste(geneListNames[iList],names(tmpDat),sep=".")
	 	} else{
	 		keepName <- tolower(names(tmpDat)) %in% tolower(c(reqID,optID))
		 	names(tmpDat)[!keepName] <- paste(geneListNames[iList],names(tmpDat)[!keepName],sep=".")
	 	}
	 	tmpDat[,paste("is",geneListNames[iList],sep=".")] <- TRUE
	  allDat <- merge(allDat,tmpDat,all=TRUE)
  }
 allDat$ZZ <- NULL
 
 # now get which list and combination of lists each ID belongs to
#  compNames <- compareNames(ccOptions)
#  geneCompVec <- .compMem(geneLists,ccOptions)
#  geneCompName <- sapply(geneCompVec,function(x) compNames[x])
 allDat <- new("mergedData", allDat, useIDName=allDatUseID)
 return(allDat)
}
	 
setMethod("mergeLists", signature=list(ccGeneList="ccGeneList",ccOptions="ccOptions"),  function(ccGeneList,ccOptions,isGene) .mergeLists(ccGeneList,ccOptions,isGene))