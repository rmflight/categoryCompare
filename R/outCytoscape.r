setMethod("ccOutCyt", signature=list(ccCompRes="ccCompareResult",ccOpts="ccOptions"),
					function(ccCompRes,ccOpts,postText="",layout=NULL,...) .ccOutCyt(ccCompRes,ccOpts,postText,layout,...))

.ccOutCyt <- function(ccCompRes,ccOpts,postText="",layout=NULL,...){
	graphName <- paste(ccCompRes@categoryName,ccCompRes@ontology,postText,sep=":")
	ccGraph <- mainGraph(ccCompRes)
	graphLayout <- ccGraph@graphData$layout
	# default is to use whatever is defined from the mainGraph, but if the user supplies a layout, then use that instead, providing an easy override of the default
  if (!(is.null(layout))){
  	graphLayout <- layout
  }
  cw <- createNetworkFromGraph(ccGraph, title = graphName, ...)
	tmpCols <- compareColors(ccOpts)
	names(tmpCols) <- NULL
 	nodeAtts <- names(nodeData(ccGraph,nodes(ccGraph)[1])[[1]])
  toolTipLoc <- grep("tooltip", nodeAtts, ignore.case=TRUE,value=FALSE)
  if (length(toolTipLoc) > 0){
  	setNodeTooltipMapping(nodeAtts[toolTipLoc[1]], network = cw)
  }
	setNodeColorMapping('fillcolor', tmpCols, tmpCols, mapping.type = "passthrough", default.color='#FF0000', network = cw)

	nodeShapes <- unique(unlist(nodeData(ccGraph,,"shape")))
	setNodeShapeMapping('shape', nodeShapes, nodeShapes, default.shape='ELLIPSE',  network = cw)
	layoutNetwork(graphLayout, network = cw)
	return(cw)
}

# used for working with categoryCompare graphs in RCytoscape
setMethod("breakEdges", signature=list(cwObject="numeric",cutoff="numeric"), function(cwObject,cutoff,edgeAtt,valDir,layout) .breakEdges(cwObject,cutoff,edgeAtt,valDir,layout))

.breakEdges <-	function(cwObject,cutoff,edgeAtt='weight',valDir='under',layout='force-directed'){
	edgeDat <- getTableColumns(table = "edge", network = cwObject)

	switch(valDir,
					under = edgeDat <- edgeDat[(as.numeric(edgeDat[,edgeAtt]) < cutoff),],
					over = edgeDat <- edgeDat[(as.numeric(edgeDat[,edgeAtt]) > cutoff),],
	)


	if (nrow(edgeDat) > 0) {
	  selectedEdges <- selectEdges(edgeDat$SUID, network = cwObject)
	  deletedEdges <- deleteSelectedEdges(network = cwObject)

	  if (!(is.null(layout)) | !(length(layout) == 0)){
	    layoutNetwork(layout, network = cwObject)
	    #layout(cwObject, layout)
	  } else {
	    layoutNetwork(network = cwObject)
	    #layout(cwObject)
	  }
	}

  message("Removed ", nrow(edgeDat), " edges from graph\n")
}


setMethod("breakEdges", signature=list(cwObject="ccCompareResult", cutoff="numeric"), function(cwObject, cutoff, edgeAtt, valDir) .breakEdgesGraph(cwObject, cutoff, edgeAtt, valDir))

.breakEdgesGraph <- function(cwObject, cutoff, edgeAtt='weight', valDir='under'){
	inGraph <- cwObject@mainGraph
	edgeW <- unlist(edgeData(inGraph,,,edgeAtt))
	switch(valDir,
				 under = delEdges <- names(edgeW)[edgeW < cutoff],
				 over = delEdges <- names(edgeW)[edgeW > cutoff],
				 )
	if (length(delEdges) > 0){
		delEdges <- strsplit(delEdges,'|', fixed=TRUE)
		fromEdges <- sapply(delEdges, function(x){x[1]})
		toEdges <- sapply(delEdges, function(x){x[2]})
		inGraph <- removeEdge(fromEdges,toEdges,inGraph)
	}
	cwObject@mainGraph <- inGraph
	message("Removed ", length(delEdges)/2, " edges from graph\n")
	return(cwObject)
}


setMethod("cwReload", signature=list(oldCW="CytoscapeWindowClass",windowName="character",ccOpts="ccOptions"), function(oldCW,windowName,ccOpts,...) .cwReload(oldCW,windowName,ccOpts,...))

# Re-connect to Cytoscape containing a graph from an old CytoscapeConnection instance
.cwReload <-	function(oldCW,windowName,ccOpts,rpcPort=9000,host="localhost"){
	newCW <- existing.CytoscapeWindow(windowName, rpcPort=rpcPort, host=host, copy.graph.from.cytoscape.to.R=FALSE)
	newCW@graph <- oldCW@graph # copy the graph from the old Cytoscape instance to the new one
	tmpCols <- compareColors(ccOpts)
	names(tmpCols) <- NULL
	setNodeColorRule(newCW, node.attribute.name='fillcolor', tmpCols, tmpCols, mode='lookup', default.color='#FF0000')
	return(newCW)
}

setMethod("resetColors", signature=list(cwObj="CytoscapeWindowClass",
																				ccOpts="ccOptions"),
					function(cwObj,ccOpts,...) .resetColors(cwObj,ccOpts,...))

.resetColors <- function(cwObj, ccOpts, node.attribute.name='fillcolor', mode='lookup'){
	tmpCols <- compareColors(ccOpts)
	names(tmpCols) <- NULL
	setNodeColorRule(cwObj, node.attribute.name=node.attribute.name, tmpCols, tmpCols, mode=mode)
}

setMethod("minNodes", signature=list(cwObj="CytoscapeWindowClass",cutoff="numeric"), function(cwObj,cutoff) .minNodes(cwObj,cutoff))

.minNodes <-	function(cwObj,cutoff){
	nodeAtts <- getAllNodeAttributes(cwObj)
	hasCount <- grep('[[:punct:]]Count',names(nodeAtts),ignore.case=TRUE)

	nCount <- length(hasCount)
	throwNode <- vector('logical',nrow(nodeAtts))
	nodeCount <- nodeAtts[,hasCount] < cutoff
	nodeCount <- apply(nodeCount,1,'sum')

	selectNodes(cwObj,names(nodeCount)[nodeCount == nCount])
	deleteSelectedNodes(cwObj)
	layoutNetwork(cwObj,'force-directed')
	#layout(cwObj, 'force-directed')
  message("Removed ", sum(nodeCount == nCount), " nodes from graph")
}

# this gets which nodes are currently selected
setMethod("cytOutNodes", signature=list(descStr="character", cwObj="CytoscapeWindowClass", saveObj="list"), function(descStr, cwObj, saveObj, outImages) .cytOutNodes(descStr, cwObj, saveObj, outImages))

setMethod("cytOutNodes", signature=list(descStr="character", cwObj="CytoscapeWindowClass", saveObj="missing"), function(descStr, cwObj, saveObj, outImages) .cytOutNodes(descStr, cwObj, saveObj, outImages))

.cytOutNodes <- function(descStr,cwObj,saveObj=vector('list',0),outImages=NULL){
	numEnt <- length(saveObj) + 1
 	currNodes <- getSelectedNodes(cwObj)
	nNodes <- length(currNodes)
 	if (nNodes == 0){
 		stop("No nodes selected!", call.=FALSE)
 	}
	if (missing(descStr) || (is.null(descStr))){
		descStr <- paste("Group",numEnt,collapse=".")
	}
 	tmpDat <- nodeData(cwObj@graph,currNodes)
  saveObj[[numEnt]] <- list(descStr=descStr,nodes=currNodes,nodeData=tmpDat)
 	if (!is.null(outImages)){
 		if (dirname(outImages) == '.'){
 			currDir <- getwd()
	 		fullPath <- file.path(currDir,outImages)
 		} else { fullPath <- outImages }
	 	dir.create(fullPath,showWarnings=FALSE)
	 	fileName <- file.path(fullPath,paste(descStr,"png",sep="."))
	 	saveImage(cwObj,fileName,"png",1)
 	}
 	return(saveObj)
}

# and then we need to get out the items annotated to those nodes (if applicable), and the data, and save it to a file if a filename is provided
setMethod("cytOutData", signature=list(saveObj='list', compareResult="ccCompareResult", mergedData="mergedData"), function(saveObj, compareResult, mergedData, orgType, fileName, displayFile) .cytOutData(saveObj, compareResult, mergedData, orgType, fileName, displayFile))

setMethod("cytOutData", signature=list(saveObj='list', compareResult="missing", mergedData="missing"), function(saveObj, compareResult=NULL, mergedData=NULL, orgType, fileName, displayFile) .cytOutData(saveObj, compareResult=NULL, mergedData=NULL, orgType, fileName, displayFile))

setMethod("cytOutData", signature=list(saveObj='list', compareResult="ccCompareResult", mergedData="missing"), function(saveObj, compareResult, mergedData=NULL, orgType, fileName, displayFile) .cytOutData(saveObj, compareResult, mergedData=NULL, orgType, fileName, displayFile))

.cytOutData <- function(saveObj,compareResult=NULL,mergedData=NULL,orgType="header",fileName=NULL,displayFile=FALSE){
	if (is.null(fileName)){
		fileName <- tempfile()
		displayFile <- TRUE
	} else {
		if (dirname(fileName) == '.'){
 			currDir <- getwd()
	 		fileName <- file.path(currDir,fileName)
 		}
	}
	if (orgType == "header"){
		outData <- .headerOutData(saveObj,compareResult,mergedData,fileName)
	} else if (orgType == "annotate"){
		outData <- .annotateOutData(saveObj,compareResult,mergedData,fileName)
	}

	if (displayFile){
		file.show(fileName,title="ccCompareResults")
	} else {
		return(outData)
	}

}

# this splits the tables up into the chunks that belong in each grouping defined by the user
.headerOutData <- function(saveObj,compareResult,mergedData,fileName){
	nSave <- length(saveObj)
 	useMerged <- TRUE # are we using a merged data table
	useAnn <- FALSE 	# are we using annotations (i.e. do we know which genes are annotated with what)
	returnDat <- vector('list', nSave)
	names(returnDat) <- sapply(saveObj, function(x){x$descStr})
 	if (is.null(mergedData)){
 		useMerged <- FALSE
 	}
 	if (is.null(compareResult)){
 		mainTable <- nodeDat2Table(saveObj[[1]]$nodeData)
		allAnn <- NULL
 		for (iSave in 2:nSave){
 			mainTable <- rbind(nodeDat2Table(saveObj[[iSave]]$nodeData))
 		}
 	} else {
 		mainTable <- compareResult@mainTable
 		allGraph <- compareResult@mainGraph
		mainTable <- addListMembership(mainTable, allGraph) # add list membership from the graph
 		allAnn <- compareResult@allAnnotation
 		useAnn <- TRUE
 	}

 	mainTable <- unique(mainTable)
	fileCon <- file(fileName,open="w+")
 	for (iSave in 1:nSave){
 		useNodes <- saveObj[[iSave]]$nodes
 		keepTable <- mainTable[match(saveObj[[iSave]]$nodes,mainTable$ID,nomatch=0),]

 		returnDat[[iSave]] <- list(AnnotationData=keepTable)

 		cat("\n\n",saveObj[[iSave]]$descStr,"\n","Annotation Data","\n",file=fileCon)
	 	write.table(keepTable,file=fileCon,sep="\t",row.names=FALSE)
		if (useAnn && useMerged) {
			tmpAnn <- allAnn[useNodes]
			tmpGenes <- unique(unlist(tmpAnn,recursive=TRUE,use.names=FALSE))
			useID <- unique(mergedData@useIDName)
			keepRow <- vector('logical',nrow(mergedData))
			for (iID in 1:length(useID)){
				keepRow <- keepRow | (mergedData[,useID[iID]] %in% tmpGenes)
			}
			keepTable <- mergedData[keepRow,]
			returnDat[[iSave]]$ItemData <- keepTable

			cat("Item Data","\n",file=fileCon)
			write.table(keepTable,file=fileCon,sep="\t",row.names=FALSE)

		}
 	}
 	close(fileCon)
 	message("Wrote file: ",fileName)
	returnDat
}

# this takes the data tables (both the annotation data and item data if available) and adds columns that indicate which user defined grouping
.annotateOutData <- function(saveObj,compareResult,mergedData,fileName){
	nSave <- length(saveObj)
 	useMerged <- TRUE
	useAnn <- FALSE
 	if (is.null(mergedData)){
 		useMerged <- FALSE
 	}
 	if (is.null(compareResult)){
 		mainTable <- nodeDat2Table(saveObj[[1]]$nodeData)
		allAnn <- NULL
 		for (iSave in 2:nSave){
 			mainTable <- rbind(nodeDat2Table(saveObj[[iSave]]$nodeData))
 		}
 	} else {
 		mainTable <- compareResult@mainTable
 		allGraph <- compareResult@mainGraph
		mainTable <- addListMembership(mainTable, allGraph)
 		allAnn <- compareResult@allAnnotation
 		useAnn <- TRUE
 	}
 	mainTable <- unique(mainTable)
	if (useAnn && useMerged){
		useID <- unique(mergedData@useIDName)
	}

 	for (iSave in 1:nSave){
 		tableName <- make.names(saveObj[[iSave]]$descStr) # create valid column names
 		mainTable[,tableName] <- FALSE
 		useNodes <- saveObj[[iSave]]$nodes
 		changeIndx <- match(saveObj[[iSave]]$nodes,mainTable$ID,nomatch=0)
 		mainTable[changeIndx,tableName] <- TRUE

		if (useAnn && useMerged){
			mergedData[,tableName] <- FALSE
			tmpAnn <- allAnn[useNodes]
			tmpGenes <- unique(unlist(tmpAnn,recursive=TRUE,use.names=FALSE))
			keepRow <- vector('logical',nrow(mergedData))
			for (iID in 1:length(useID)){
				keepRow <- keepRow | (mergedData[,useID[iID]] %in% tmpGenes)
			}
			mergedData[keepRow,tableName] <- TRUE
		}
 	}

	fileCon <- file(fileName,open="w+")
	cat("Annotation Table","\n",file=fileCon)
	write.table(mainTable,file=fileCon,row.names=FALSE)
	if (useAnn && useMerged){
		cat("\n","Item Table","\n",file=fileCon)
		write.table(mergedData,file=fileCon,row.names=FALSE)
	}
	close(fileCon)
 	message("Wrote file: ",fileName)

	returnDat <- list(AnnotationData=mainTable, ItemData=mergedData)

}

nodeDat2Table <- function(nodeDat){
	col.Names <- c("ID", names(nodeDat[[1]]))
 	row.Names <- names(nodeDat)
 	tmpDat <- matrix(0,length(row.Names),length(col.Names))
 	tmpDat <- as.data.frame(tmpDat, stringsAsFactors=FALSE)
 	for (iRow in 1:length(row.Names)){
 		tmpDat[iRow,] <- c(row.Names[iRow],unlist(nodeDat[[iRow]]))
 	}
 	names(tmpDat) <- col.Names
 	return(tmpDat)
}

# this function simply adds the "listMembership" to the mainTable
addListMembership <- function(mainTable, allGraph){
	allNodes <- nodes(allGraph)
	listMem <- unlist(nodeData(allGraph, allNodes, "listMembership"))
	tableID <- mainTable$ID

	matchID2Node <- match(allNodes, tableID, nomatch=0)

	mainTable$listMembership <- "NA"
	mainTable$listMembership[matchID2Node] <- listMem
	return(mainTable)
}
