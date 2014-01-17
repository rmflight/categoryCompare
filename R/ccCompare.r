# ccCompare - here we have the functions that actually generate comparisons between the enriched annotations determined
# for each gene list. operates on a ccEnrichCollection object, with a ccOptions object, and a directory name to save everything to.
# Must be able to dispatch different methods depending on whether or not the sub-objects are GO or KEGG or generic.

# Note: need to make this a generic eventually and then use a method definition so that we only operate on the correct classes

setMethod("ccCompare", signature=list(ccEnrichResult="ccEnrichCollection",ccOptions="ccOptions"),
					function(ccEnrichResult, ccOptions) .ccCompare(ccEnrichResult, ccOptions))

.ccCompare <- function(ccEnrichResult, ccOptions){
	
	# apply ccCompare to each sub-object, using the 
  allCompare <- lapply(ccEnrichResult, ccCompare, ccOptions)
  allCompare <- new("ccCompareCollection", allCompare)
}

setMethod("ccCompare", signature=list(ccEnrichResult="GOccEnrichResult", ccOptions="ccOptions"),
					function(ccEnrichResult, ccOptions) .ccCompareGO(ccEnrichResult, ccOptions))

.ccCompareGO <- function(ccEnrichResult, ccOptions){
  pieData <- "NA"
  annOpt <- annStatus() # what can we do. This will possibly change the type of results we can generate
	if (!annOpt$godb){
		stop('GO.db needs to be loaded for GOccEnrichResult objects!')
	}
  GOstrings <- c(BP = GOBPPARENTS, MF=GOMFPARENTS, CC=GOCCPARENTS)
  GOnames <- names(GOstrings)
  lists <- names(ccEnrichResult)
  
  useLists <- listNames(ccOptions)
  
  if (sum(lists %in% useLists) == 0){
  	stop('listNames defined in ccOptions do not match any of the names in the GOccEnrichResult object!')
  } else {
  	ccEnrichResult <- ccEnrichResult[useLists]
  	lists <- names(ccEnrichResult)
  	nList <- length(lists)
  }
  
  ccGOType <- ccEnrichResult@ontology
  
  # Extract the results from the ccEnrichResult object
  extRes <- .extractRes(ccEnrichResult)
  allTable <- extRes$allTable
  allNodes <- extRes$allNodes
  allRes <- extRes$allRes
  sigID <- extRes$sigID
	
  rm(extRes) # free up some memory
   
  # get a more descriptive item for each ID
  idDat <- mget(allTable$ID, envir=GOTERM, ifnotfound=NA)
  idDat <- lapply(idDat, Term)
  idDat <- unlist(idDat, use.names=FALSE)
  allTable$Desc <- idDat
  
  # now rearrange a few things in the table
  oldIndx <- match(lists,names(allTable))
  newIndx <- seq(2,2+length(oldIndx)-1,1)
  allTable <- .moveTable(allTable,oldIndx,newIndx)
  allTable <- .moveTable(allTable,match("Desc",names(allTable)),2)
  
  # get the genes annotated to each term by list membership
  allAnn <- .annGenesComp(allRes,ccOptions)
  
  # now lets generate the graph based on the significant nodes
  tmpList <- vector('list',nList)
  nodeListMem <- vector('list',0)
  if (graphType(ccEnrichResult) %in% "hierarchical"){
    useRes <- GOstrings[[GOnames[names(GOstrings) %in% ccGOType]]]
    allNodes <- vector('character', 0)
    for (iList in 1:nList){
      allRes[[iList]]$indGO <- nodes(GOGraph(allRes[[iList]]$sigIDs, useRes))
      allNodes <- c(allNodes,allRes[[iList]]$indGO)
      nodeListMem[[iList]] <- allRes[[iList]]$indGO
    }
    allNodes <- unique(allNodes)
    allGraph <- GOGraph(allNodes, useRes)
		allNodes <- nodes(allGraph)
		names(allNodes) <- NULL
		nodes(allGraph) <- allNodes
    allGraph@graphData$layout <- "hierarchical"
  } else if (graphType(ccEnrichResult) %in% "overlap") {
  	
    for (iList in 1:nList){
      nodeListMem[[iList]] <- allRes[[iList]]$sigIDs
      allRes[[iList]]$indGO <- allRes[[iList]]$sigIDs
    }
    
		allAnnotation <- unique(sapply(allRes, function(x){x$annotation}))
		
		nAnn <- length(allAnnotation)
		nSig <- length(sigID)
		nodeGeneMap <- vector("list",nSig) 
		names(nodeGeneMap) <- sigID
    if (nAnn > 1){
			for (iAnn in 1:nAnn){
				tmpAnn <- getGO2ALLEGS(sigID,allAnnotation[iAnn])
	 			for (iSig in 1:nSig){
	 				nodeGeneMap[[iSig]] <- unique(c(nodeGeneMap[[iSig]], tmpAnn[[iSig]]))
	 			}	 		
			}
    } else {
    	nodeGeneMap <- getGO2ALLEGS(sigID, allAnnotation)
    }
    
		
    allGraph <- createGraph2(sigID,nodeGeneMap,ccEnrichResult@overlapType)
    allGraph@graphData$layout <- "force-directed"
    
  }
  names(nodeListMem) <- lists
    
  # filter out the 'all' go node (depending on how we generated the graph it might be there
  allNodes <- nodes(allGraph)
  allNodes <- allNodes[!(allNodes %in% 'all')]
  nNodes <- length(allNodes)
  
  allGraph <- .initGraphAtts(allGraph,allTable)
  
  idDat <- mget(allNodes, envir=GOTERM, ifnotfound=NA)
  idDat <- lapply(idDat, Term)
  idDat <- unlist(idDat, use.names=FALSE)
  nodeData(allGraph, allNodes, attr="Desc") <- idDat
  
  if (colorType(ccOptions) == "solid"){
    nodeCompVec <- .compMem(nodeListMem,ccOptions)
    nodeCompVec <- nodeCompVec[match(allNodes,names(nodeCompVec),nomatch=0)]
    
    # now that we know which of the lists we belong to, we can set up some attributes
    
    nodeData(allGraph, allNodes, attr="fillcolor") <- sapply(nodeCompVec, function(x){compareColors(ccOptions)[x]}) # this is why we are supposed to do the induced graph from each, and then combine them.
    nodeData(allGraph, allNodes, attr="listMembership") <- sapply(nodeCompVec, function(x){compareNames(ccOptions)[x]})
    
    # check fillColor and listMembership, if any are missing, set them to NA
    tmpMember <- sapply(nodeData(allGraph,allNodes,"listMembership"),length)
    nodeData(allGraph,names(tmpMember)[tmpMember == 0],"listMembership") <- 'NA'
    nodeData(allGraph,names(tmpMember)[tmpMember == 0],"fillcolor") <- 'NA'
    
    nodeData(allGraph, allNodes, attr="compIndx") <- nodeCompVec # which comparison are we (if we need to access that again)
    nodeData(allGraph, allNodes[allNodes %in% sigID], attr="isSig") <- as.character(TRUE)
    nodeData(allGraph, allNodes, attr="toolTip") <- paste(unlist(nodeData(allGraph, allNodes, attr="listMembership")), allNodes,unlist(nodeData(allGraph, allNodes, attr="Desc")), sep=" <br> ")
  } else if (colorType(ccOptions) == "pie"){
    pieData <- .genPieMatrix(allNodes, allRes, ccOptions)
  }
   # and reorder to be in the same order as the nodes in the graph
  
  
    
  # only do this if we are looking at just the overlap between lists without GO context
#   if (allGraph@graphData$layout == "neato") {
#     nodeData(allGraph, allNodes[allNodes %in% lists], attr="shape") <- "box"
#   } else {
#     nodeData(allGraph, allNodes[allNodes %in% sigID], attr="shape") <- "box"
#   }
  # All of these attributes are stored now in such a way that the user or one of 
  # our functions can access them later, even if Rgraphviz is not being used.
  # maybe can use RCytoscape to view this stuff?
  
  
  
  # now we have a graph, and a table, and the annotated genes to each entry in the table. So lets give all that back
  # to the user. 
	returnData <- new("ccCompareResult", mainGraph=allGraph, mainTable=allTable, allAnnotation=allAnn, categoryName="GO", ontology=ccGOType, pieData=pieData)
  #returnData <- list(graphs=list(mainGraph=allGraph), mainTable=allTable, allAnnotation=allAnn)
  return(returnData)
}

setMethod("ccCompare", signature=list(ccEnrichResult="KEGGccEnrichResult",ccOptions="ccOptions"),
					function(ccEnrichResult, ccOptions) .ccCompareKEGG(ccEnrichResult,ccOptions))
.ccCompareKEGG <- function(ccEnrichResult, ccOptions){
  pieData <- "NA"
  annOpt <- annStatus() # what can we do. This will possibly change the type of results we can generate
	if (!annOpt$keggdb){
		stop("KEGG.db must be loaded to continue!")
	}
  summary <- getGeneric("summary")
  lists <- names(ccEnrichResult)
  
  useLists <- listNames(ccOptions)
  
  if (sum(lists %in% useLists) == 0){
  	stop('listNames defined in ccOptions do not match any of the names in the KEGGccEnrichResult object!')
  } else {
  	ccEnrichResult <- ccEnrichResult[useLists]
  	lists <- names(ccEnrichResult)
  	nList <- length(lists)
  }
  
	extRes <- .extractRes(ccEnrichResult)
  allTable <- extRes$allTable
  allNodes <- extRes$allNodes
  allRes <- extRes$allRes
  sigID <- extRes$sigID
  rm(extRes) # free up some memory
   
  # get a more descriptive item for each ID
  idDat <- unlist(mget(allTable$ID, envir=KEGGPATHID2NAME, ifnotfound=NA))
  allTable$Desc <- idDat
  
  # now rearrange a few things in the table
  oldIndx <- match(lists,names(allTable))
  newIndx <- seq(2,2+length(oldIndx)-1,1)
  allTable <- .moveTable(allTable,oldIndx,newIndx)
  allTable <- .moveTable(allTable,match("Desc",names(allTable)),2)
  
  # get the genes annotated to each term by list membership
  allAnn <- .annGenesComp(allRes,ccOptions)
  
  # now lets generate the graph based on the significant nodes
	
	nodeListMem <- vector('list',nList)
  for (iList in 1:nList){
    nodeListMem[[iList]] <- allRes[[iList]]$sigIDs
  }
  names(nodeListMem) <- lists
		
	# this allows us to have two different organisms and compare them directly
	allAnnotation <- unique(sapply(allRes, function(x){x$annotation}))
	nAnn <- length(allAnnotation)
	nSig <- length(sigID)
	nodeGeneMap <- vector("list",nSig) 
	names(nodeGeneMap) <- sigID
	for (iAnn in 1:nAnn){
		tmpAnn <- getPATH2EG(sigID,allAnnotation[iAnn])
 		for (iSig in 1:nSig){
 			nodeGeneMap[[iSig]] <- unique(c(nodeGeneMap[[iSig]], tmpAnn[[iSig]]))
 		}	 		
	}
	
  allGraph <- createGraph2(sigID,nodeGeneMap,ccEnrichResult@overlapType)
  allGraph@graphData$layout <- "force-directed"
	
  # Note that we replace allNodes now based on what is now in the graph  
  allNodes <- nodes(allGraph)
  allNodes <- allNodes[!(allNodes %in% 'all')]
  nNodes <- length(allNodes)
  
	nodeCompVec <- .compMem(nodeListMem,ccOptions)
  nodeCompVec <- nodeCompVec[match(allNodes,names(nodeCompVec),nomatch=0)] # and reorder to be in the same order as the nodes in the graph
  
  # now that we know which of the lists we belong to, we can set up some attributes
  allGraph <- .initGraphAtts(allGraph,allTable)

	idDat <- unlist(mget(allNodes, envir=KEGGPATHID2NAME, ifnotfound=NA))
	nodeData(allGraph, allNodes, attr="Desc") <- idDat
	
  if (colorType(ccOptions) == "solid"){
    nodeData(allGraph, allNodes, attr="fillcolor") <- sapply(nodeCompVec, function(x){compareColors(ccOptions)[x]}) # this is why we are supposed to do the induced graph from each, and then combine them.
    nodeData(allGraph, allNodes, attr="listMembership") <- sapply(nodeCompVec, function(x){compareNames(ccOptions)[x]})
    
    # check fillColor and listMembership, if any are missing, set them to NA
    tmpMember <- sapply(nodeData(allGraph,allNodes,"listMembership"),length)
    nodeData(allGraph,names(tmpMember)[tmpMember == 0],"listMembership") <- 'NA'
    nodeData(allGraph,names(tmpMember)[tmpMember == 0],"fillcolor") <- 'NA'
    
    nodeData(allGraph, allNodes, attr="compIndx") <- nodeCompVec # which comparison are we (if we need to access that again)
    nodeData(allGraph, allNodes[allNodes %in% sigID], attr="isSig") <- as.character(TRUE)
    nodeData(allGraph, allNodes, attr="toolTip") <- paste(unlist(nodeData(allGraph, allNodes, attr="listMembership")), allNodes,unlist(nodeData(allGraph, allNodes, attr="Desc")), sep=" <br> ") 
  } else if (colorType(ccOptions) == "pie"){
    pieData <- .genPieMatrix(allNodes, allRes, ccOptions)
  }
  
  
    
  # only do this if we are looking at just the overlap between lists without GO context
#   if (allGraph@graphData$layout == "neato") {
#     nodeData(allGraph, allNodes[allNodes %in% lists], attr="shape") <- "box"
#   } else {
#     nodeData(allGraph, allNodes[allNodes %in% sigID], attr="shape") <- "box"
#   }
  # All of these attributes are stored now in such a way that the user or one of 
  # our functions can access them later, even if Rgraphviz is not being used.
  # maybe can use RCytoscape to view this stuff?
  
  
  
  # now we have a graph, and a table, and the annotated genes to each entry in the table. So lets give all that back
  # to the user. 
  returnData <- new("ccCompareResult", mainGraph=allGraph, mainTable=allTable, allAnnotation=allAnn, categoryName="KEGG", pieData=pieData)
	#returnData <- list(graphs=list(mainGraph=allGraph), mainTable=allTable, allAnnotation=allAnn)
  return(returnData)
}

setMethod("ccCompare", signature=list(ccEnrichResult="GENccEnrichResult", ccOptions="ccOptions"),
					function(ccEnrichResult, ccOptions) .ccCompareGeneric(ccEnrichResult, ccOptions))
.ccCompareGeneric <- function(ccEnrichResult, ccOptions){
  ### we need to add a slot, overlapType, to this and the others that sets the overlap to use based on the categoryName, and allows the user to override it if they want
	annOpt <- annStatus() # what can we do. This will possibly change the type of results we can generate
	annStr <- NULL
  graphT <- graphType(ccEnrichResult)
  hasMap <- FALSE
  if (length(geneAnnMapping(ccEnrichResult)) > 0){
  	hasMap <- TRUE
  }
  if (!hasMap){
  	switch(categoryName(ccEnrichResult),
  				 KEGG = annStr <- "PATH2EG",
  				 GO = annStr <- "GO2ALLEG",
  				 NULL) # put in null as an option, should allow simple control of behavior of other functions
  	if ((!annOpt$keggdb) && (categoryName(ccEnrichResult) == "KEGG")){
  		stop("KEGG requested, but KEGG.db not loaded and no geneAnnMapping supplied!")	 
  	} else if ((!annOpt$godb) && (categoryName(ccEnrichResult) == "GO")){
  		stop("GO requested, but GO.db not loaded and no geneAnnMapping supplied!")
  	}
  } 
	if ((graphT == "hierarchical") && (categoryName(ccEnrichResult) != "GO")){
		stop("Hierarchical GO layout only supported for GO!")
	}
	
  lists <- names(ccEnrichResult)
  
  useLists <- listNames(ccOptions)
  
  if (sum(lists %in% useLists) == 0){
  	stop('listNames defined in ccOptions do not match any of the names in the ccEnrichResult object!')
  } else {
  	ccEnrichResult <- ccEnrichResult[useLists]
  	lists <- names(ccEnrichResult)
  	nList <- length(lists)
  }
 
  sigIDs <- unique(unlist(sapply(ccEnrichResult,function(x){sigID(x)})))
  if (!hasMap){
  	tmpDesc <- .getDesc(sigIDs,categoryName(ccEnrichResult))
  } else {
    tmpDesc <- rep("", length(sigIDs))
    names(tmpDesc) <- sigIDs
  }
	allTable <- data.frame(ID=sigIDs,Desc=tmpDesc)
  
  # now lets generate the graph based on the significant nodes
	
	nodeListMem <- vector('list',nList)
  for (iList in 1:nList){
    nodeListMem[[iList]] <- sigID(ccEnrichResult[[iList]])
  }
  names(nodeListMem) <- lists
 
  # this allows us to have two different organisms and compare them directly
  if (!hasMap){
  	allAnnotation <- unique(sapply(ccEnrichResult, function(x){annotation(x)}))
		nAnn <- length(allAnnotation)
		nSig <- length(sigIDs)
		nodeGeneMap <- vector("list",nSig) 
		names(nodeGeneMap) <- sigIDs
		for (iAnn in 1:nAnn){
			tmpAnn <- getAnnotation(sigIDs,allAnnotation[iAnn],annStr, F)
	 		for (iSig in 1:nSig){
	 			nodeGeneMap[[iSig]] <- unique(c(nodeGeneMap[[iSig]], tmpAnn[[iSig]]))
	 		}	 		
		}
  } else {
   	nodeGeneMap <- geneAnnMapping(ccEnrichResult)[sigIDs]
  }
	
	 #browser(expr=TRUE)
	if (graphT == "overlap"){
		allGraph <- createGraph2(sigIDs,nodeGeneMap,ccEnrichResult@overlapType)
		allGraph@graphData$layout <- "force-directed"
		nodeCompVec <- .compMem(nodeListMem,ccOptions)
	} else if (graphT == "membership"){
		allGraph <- createGraph(nodeListMem)
		allGraph@graphData$layout <- "attribute-circle"
		nodeCompVec <- .compMem(nodeListMem,ccOptions)
		tmpVec <- unlist(ccOptions@compareIndx[names(ccOptions@compareIndx) %in% listNames(ccOptions)])
		nodeCompVec <- c(nodeCompVec,tmpVec)
	}
  
	allGraph <- .initGraphAtts(allGraph,allTable)
	allNodes <- nodes(allGraph)
	nodeCompVec <- nodeCompVec[match(allNodes,names(nodeCompVec),nomatch=0)]
	
	# idDat <- allTable$Desc[match(allTable$ID, allNodes, nomatch=0)]
	nodeData(allGraph, allNodes, attr="Desc") <- tmpDesc[allNodes]

  nodeData(allGraph, allNodes, attr="fillcolor") <- sapply(nodeCompVec, function(x){compareColors(ccOptions)[x]}) # this is why we are supposed to do the induced graph from each, and then combine them.
  nodeData(allGraph, allNodes, attr="listMembership") <- sapply(nodeCompVec, function(x){compareNames(ccOptions)[x]})
  
  # check fillColor and listMembership, if any are missing, set them to NA
  tmpMember <- sapply(nodeData(allGraph,allNodes,"listMembership"),length)
  nodeData(allGraph,names(tmpMember)[tmpMember == 0],"listMembership") <- 'NA'
  nodeData(allGraph,names(tmpMember)[tmpMember == 0],"fillcolor") <- 'NA'
  
  nodeData(allGraph, allNodes, attr="compIndx") <- nodeCompVec # which comparison are we (if we need to access that again)
  nodeData(allGraph, allNodes[allNodes %in% sigIDs], attr="isSig") <- as.character(TRUE)
	nodeData(allGraph, allNodes, attr="toolTip") <- paste(unlist(nodeData(allGraph, allNodes, attr="listMembership")), allNodes,unlist(nodeData(allGraph, allNodes, attr="Desc")), sep=" <br> ")
	
	# only do this if we are looking at just the overlap between lists without considering gene overlap
  if (graphT == "membership") {
    nodeData(allGraph, allNodes[allNodes %in% lists], attr="shape") <- "trapezoid"
  } 
  
  returnData <- new("ccCompareResult", mainGraph=allGraph, mainTable=allTable, categoryName=categoryName(ccEnrichResult))
	#returnData <- list(graphs=list(mainGraph=allGraph), mainTable=allTable, allAnnotation=allAnn)
  return(returnData)
  	
}

.getDesc <- function(sigIDs,catType){
	if (catType == "KEGG"){
		tmpDesc <- getAnnotation(sigIDs,"KEGG.db","PATHID2NAME",TRUE)
	} else if (catType == "GO") {
		tmpDesc <- getAnnotation(sigIDs,"GO.db", "GOTERM", FALSE)
		tmpDesc <- sapply(tmpDesc,function(x) {
			if(length(x) > 0){
				Term(x)
			} else {
				NA
			}
		})
	} else { tmpDesc <- vector("character", length(sigIDs))}
	return(tmpDesc)
}  	
  	
.moveTable <- function(dataFrame,oldIndx,newIndx){
  nCol <- ncol(dataFrame)
  allIndx <- seq(1,nCol,1)
  
  nChng <- length(oldIndx)
  for (iChng in 1:nChng){
    # chop off from newIndx down
    oldS <- allIndx[newIndx[iChng]:nCol]
    oldS <- oldS[!(oldS == oldIndx[iChng])]
    newS <- allIndx[1:newIndx[iChng]-1]
    allIndx <- c(newS,oldIndx[iChng],oldS)
  }
  dataFrame <- dataFrame[,allIndx]
  return(dataFrame)
  
}


# various types of similarity coefficients to calculate for the nodes.

jaccardCoef <- function(n1, n2){
  (length(intersect(n1,n2))) / (length(union(n1,n2)))
}

overlapCoef <- function(n1, n2){
  (length(intersect(n1,n2))) / (min(c(length(n1),length(n2))))
}

combinedCoef <- function(n1, n2){
  oCoef <- overlapCoef(n1, n2)
  jCoef <- jaccardCoef(n1, n2)
  return((0.5*oCeof) + (0.5*jCoef))
}


createGraph2 <- function(nodeList,nodeGeneMap,overlapType){
  # graph creation is based on:
  # Merico D, Isserlin R, Stueker O, Emili A, Bader GD, 2010 
  # Enrichment Map: A Network-Based Method for Gene-Set Enrichment Visualization and Interpretation. PLoS ONE 5(11): e13984. doi:10.1371/journal.pone.0013984
  
  # stop code at the spot where errors may likely creep in.
  if (length(nodeList) == 0){
  	warning("Created a graph with zero nodes and zero edges!", call.=F)
  	return(new("graphNEL", nodes=character(0), edgemode="directed"))
  }
  
  if (length(nodeList) == 1){
  	warning("Created a graph with one node and zero edges!", call.=F)
  	return(new("graphNEL", nodes=nodeList, edgemode="directed"))
  }
  
  nGenesNode <- sapply(nodeGeneMap,'length')
  keepNodes <- (nGenesNode >= 10) & (nGenesNode <= 500)
  nodeList <- nodeList[keepNodes]
  nodeGeneMap <- nodeGeneMap[keepNodes]
  nNodes <- length(nodeList)

  graphDat <- new("graphNEL", nodes=nodeList, edgemode="directed") # this is because we will be visualizing in Cytoscape, and really, we don't actually need the dual information.
  
  allComp <- expand.grid(seq(1,nNodes),seq(1,nNodes))
  allComp <- allComp[(allComp[,2] > allComp[,1]),] # keep only where second is greater than first
  allComp <- as.matrix(allComp)
  
  # now go through allComp, using them as indices into the matrix
  corAll <- sapply(seq(1,nrow(allComp)), function(x){
  	doComp <- allComp[x,]
  	n1 <- nodeGeneMap[[doComp[1]]]
  	n2 <- nodeGeneMap[[doComp[2]]]
  	
  	useC <- switch(overlapType,
                   overlap=overlapCoef(n1,n2),
                   jaccard=jaccardCoef(n1,n2),
                   combined=combinedCoef(n1,n2))
  	useC
  })
  
  # get rid of anything that was completely zero
  notZero <- corAll != 0
  allComp <- allComp[notZero,]
  corAll <- corAll[notZero]
  fromEdge <- nodeList[allComp[,1]]
  toEdge <- nodeList[allComp[,2]]
  
  graphDat <- addEdge(fromEdge, toEdge, graphDat, corAll)
  
}


createGraph <- function(nodeList){

  nMain <- length(nodeList)
	mainNodes <- names(nodeList)	
	nodeDat <- names(nodeList)
	nodeDat <- c(nodeDat, unique(unlist(nodeList), recursive=TRUE, use.names=FALSE))
	nNode <- length(nodeDat)

	baseAM <- matrix(data=0, nrow=nNode, ncol=nNode)
	rownames(baseAM) <- colnames(baseAM) <- nodeDat

	for (iMain in 1:nMain){
		nSub <- length(nodeList[[iMain]])
		mainNam <- mainNodes[iMain]
		for (iSub in 1:nSub){
			subName <- nodeList[[iMain]][iSub]
			baseAM[mainNam,subName] <- 1
			baseAM[subName,mainNam] <- 1
		}
	}

	graphDat <- new('graphAM', adjMat=baseAM)
	graphDat

}

.compMem <- function(memList,ccOptions){
  allItem <- unique(unlist(memList,use.names=FALSE))
  allItem <- allItem[!(allItem %in% 'all')] # in case we have the GO term "all" in our list, we don't want it
  nItem <- length(allItem)
  allList <- names(memList)
  matchList <- match(allList,listNames(ccOptions)) # want this so we grab the lists in the same order as was defined in ccOptions. They should be in the same order, but just in case
  nList <- length(memList)
  isInList <- matrix(data=NA, nrow=length(allItem), ncol=nList)
  for (iList in 1:nList){
    isInList[,iList] <- allItem %in% memList[[matchList[iList]]]
  }
  
  compIndx <- compareIndx(ccOptions)
  nComp <- length(compIndx)
  compVec <- vector('integer',nItem) 
    
  for (iComp in 1:nComp){
    truVec <- vector('logical', nList)
    truVec[compIndx[[iComp]]] <- TRUE
    
    matchVec <- apply(isInList,1,function(x){truVec == x})
    matchVec <- t(matchVec)
    matchVec <- apply(matchVec,1,sum)
    matchVec <- matchVec == nList
    compVec[matchVec] <- iComp
  }
  names(compVec) <- allItem
  return(compVec)
}

.annGenesComp <- function(allRes,ccOptions){
  lists <- names(allRes)
  nList <- length(lists)
  geneLists <- vector('list',nList)
  for (iList in 1:nList){
    geneLists[[iList]] <- allRes[[iList]]$genes
  }
  names(geneLists) <- lists
  geneCompVec <- .compMem(geneLists,ccOptions) # good to here
  
  # now replace the annGenes list by a list that includes which of the comparisons the genes belong to
  compNames <- compareNames(ccOptions)
  geneCompName <- sapply(geneCompVec,function(x){compNames[x]})  
  names(geneCompName) <- names(geneCompVec)
  nComp <- length(compNames)
  
  newAnnGenes <- vector('list',0)
  for (iList in 1:nList){
    oldNames <- names(newAnnGenes)
    currNames <- names(allRes[[iList]]$annGenes)
    matchNames <- match(currNames, oldNames)
    nCurr <- length(currNames)
    
    for (iCurr in 1:nCurr){
      currGenes <- allRes[[iList]]$annGenes[[iCurr]]
      tmpGenes <- vector('list',nComp)
      names(tmpGenes) <- compNames
      if(is.na(matchNames[iCurr])){
        
        for (iComp in 1:nComp){
          posGenes <- names(geneCompVec[geneCompName == compNames[iComp]])
          tmpGenes[[iComp]] <- intersect(posGenes,currGenes) # genes annotated to that term in this comparison
        }
        nAnn <- length(newAnnGenes)
        newAnnGenes[[nAnn+1]] <- tmpGenes
        names(newAnnGenes)[nAnn+1] <- currNames[iCurr]
        
      } else {
        oldGenes <- newAnnGenes[[matchNames[iCurr]]]
        for (iComp in 1:nComp){
          posGenes <- names(geneCompName[geneCompName == compNames[iComp]])
          tmpGenes[[iComp]] <- unique(c(oldGenes[[iComp]],intersect(posGenes,currGenes)))
        }
        newAnnGenes[[matchNames[iCurr]]] <- tmpGenes
      }
      
    }
  }
  return(newAnnGenes)
}

.extractRes <- function(ccEnrichResult){
	summary <- getGeneric("summary")
	lists <- names(ccEnrichResult)
	nList <- length(ccEnrichResult)
	
	sigID <- vector('character',0)
	allRes <- vector('list', nList)
	allNodes <- vector('character',0)
	for (iList in 1:nList){
    tmpRes <- vector('list', 4)
    tmpSum <- summary(ccEnrichResult[[iList]])
    tmpRes[[1]] <- tmpSum$ID # just the significant IDs
    tmpRes[[2]] <- summary(ccEnrichResult[[iList]], pvalue=1, pType='pval', minCount=0) # everything in the table
    tmpRes[[3]] <- geneIdsByCategory(ccEnrichResult[[iList]]) # which genes from the list are annotated to each term
    tmpRes[[4]] <- geneIds(ccEnrichResult[[iList]]) # the genes that were input to the enrichment calculation
		tmpRes[[5]] <- annotation(ccEnrichResult[[iList]]) # what annotation was used to generate the data
    
    names(tmpRes) <- c('sigIDs','allData','annGenes','genes','annotation')
    allRes[[iList]] <- tmpRes
    sigID <- c(sigID,tmpRes[[1]])
  
    allNodes <- c(allNodes, tmpRes[[1]])
  }
  names(allRes) <- lists
  allNodes <- unique(allNodes)
  sigID <- unique(sigID)
  
  tmpTable <- allRes[[1]]$allData
  nCol <- ncol(tmpTable)
  names(tmpTable)[2:nCol] <- paste(lists[1],names(tmpTable)[2:nCol],sep='.')
	
	if (nrow(tmpTable) > 0){
  	tmpTable[lists[1]] <- TRUE # this needs better behavior when there is no results for that list
	}
	
	# what if we remove those that don't have any significant results first
  
  allTable <- tmpTable
  
  if (nList > 1){
		for (iList in 2:nList){
			tmpTable <- allRes[[iList]]$allData
			nCol <- ncol(tmpTable)
			names(tmpTable)[2:nCol] <- paste(lists[iList],names(tmpTable)[2:nCol],sep='.')
			if (nrow(tmpTable) > 0) { tmpTable[lists[iList]] <- TRUE } # only put it in if there is more than one row, otherwise it will die
			allTable <- merge(allTable, tmpTable, by.x="ID", by.y="ID", all=TRUE)
		}
  }
  
  return(list(allRes=allRes,allNodes=allNodes,allTable=allTable,sigID=sigID))
}

# add attributes to the graph and add on the table data to the graph as well
.initGraphAtts <- function(allGraph,tableDat){
	typeConv <- c('STRING','INTEGER','DOUBLE','STRING')
  names(typeConv) <- c('character','integer','numeric','logical')
  
	# Note that we add an attribute so that we are RCytoscape friendly for almost no extra work
  nodeDataDefaults(allGraph, "shape") <- "ellipse"
  attr(nodeDataDefaults(allGraph, "shape"), "class") <- "STRING"
  nodeDataDefaults(allGraph, "Desc") <- ""
	attr(nodeDataDefaults(allGraph, "Desc"), "class") <- "STRING"
  nodeDataDefaults(allGraph, "listMembership") <- ""
  attr(nodeDataDefaults(allGraph, "listMembership"), "class") <- "STRING"
  nodeDataDefaults(allGraph, "compIndx") <- ""
  attr(nodeDataDefaults(allGraph, "compIndx"), "class") <- "STRING"
  nodeDataDefaults(allGraph, "fillcolor") <- ""
  attr(nodeDataDefaults(allGraph, "fillcolor"), "class") <- "STRING"
  nodeDataDefaults(allGraph, "toolTip") <- ""
  attr(nodeDataDefaults(allGraph, "toolTip"), "class") <- "STRING"
  nodeDataDefaults(allGraph, "isSig") <- "FALSE"
	attr(nodeDataDefaults(allGraph, "isSig"), "class") <- "STRING"
  
	tmpDats <- names(edgeDataDefaults(allGraph))
  if (!("weight" %in% tmpDats)){
  	edgeDataDefaults(allGraph,"weight") <- 1
  }
  attr(edgeDataDefaults(allGraph, "weight"), "class") <- "DOUBLE"
  
  allNodes <- nodes(allGraph)
  # rearrange the table to have the same ordering as the nodes in the graph
  tableDat <- tableDat[match(allNodes,tableDat$ID,nomatch=0),]
  allNodes <- allNodes[match(tableDat$ID,allNodes,nomatch=0)]
  # now work through the columns of the table, excluding the ID column
  excludeCols <- c('ID','Desc')
  keepCols <- names(tableDat)[!(names(tableDat) %in% excludeCols)]
  nKeep <- length(keepCols)
  if (nKeep > 0){
	  for (iKeep in 1:nKeep){
	  	tmpDat <- tableDat[,keepCols[iKeep]]
			tmpDat[is.na(tmpDat)] <- -1
			tmpDat[is.infinite(tmpDat)] <- -1
	  	colClass <- class(tmpDat)
	  	datName <- keepCols[iKeep]
	  	
	  	if (colClass == 'logical'){
	  		nodeDataDefaults(allGraph, datName) <- "FALSE"
	  		attr(nodeDataDefaults(allGraph, datName), "class") <- "STRING"
	  		nodeData(allGraph,allNodes,datName) <- as.character(tmpDat)
	  	} else {
	  		nodeDataDefaults(allGraph, datName) <- as(-1, colClass)
			  tmpType <- typeConv[colClass]; names(tmpType) <- NULL;
	 			attr(nodeDataDefaults(allGraph, datName), "class") <- typeConv[colClass]
				nodeData(allGraph,allNodes,datName) <- tmpDat
	  	}
	  }
  }
  return(allGraph)
}

.genPieMatrix <- function(allNodes, allRes, ccOptions){
  nNode <- length(allNodes)
  pieMatrix <- do.call(cbind, lapply(names(allRes), function(inName){
    baseCol <- rep(ccOptions@unsaturatedColor[inName], nNode)
    changeIndex <- which(allNodes %in% allRes[[inName]]$sigID)
    baseCol[changeIndex] <- compareColors(ccOptions)[inName]
    baseCol
  }))
  rownames(pieMatrix) <- allNodes
  pieDir <- file.path(getwd(), "pieGraphs")
  dir.create(pieDir)
  
  pieArea <- rep(1 / length(names(allRes)), length(names(allRes)))
  names(pieArea) <- ""
  pieFiles <- sapply(rownames(pieMatrix), function(inRow){
    filename <- file.path(pieDir, paste(sub(":", "_", inRow, fixed=TRUE), ".png", sep=""))
    png(filename=filename, bg="transparent")
    pie(pieArea, col=pieMatrix[inRow,], clockwise=TRUE)
    dev.off()
    return(filename)
  })
  return(pieFiles)
}

#  This may be supported later, but not right now. Have to think about how to support multiple organisms
#  May just refuse to do it on multiple organisms
# keggSub <- function(ccResults,allRes,allGraph,keggData){
#   lists <- names(ccResults)
#   nList <- length(lists)
#   
#   for (iList in 1:nList){
#     # check if the organism is supported
#     
#   
#   
# }