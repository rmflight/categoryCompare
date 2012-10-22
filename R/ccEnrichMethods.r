# Here we define the methods for a ccGeneList object that actually generate the enrichments and spit out the 
# enriched annotations. And then we will save this to an enrichment object that we will pass along with
# ccOptions to generate the actual output.

# this is the overarching function that controls which sub-functions are called
# Note: should be made into a generic and then a method defined.
setMethod("ccEnrich", "ccGeneList", 
	function(ccGeneList){
	  goType <- c("BP","MF","CC")
	  keggType <- "KEGG"
	  anyType <- "ANY"
	  annOpt <- annStatus()
	  
	  ccNames <- names(ccGeneList)
	  
	  # some checks that we can do what was requested. This should eventually work for custom annotations, but not sure how that will work yet.
	  testCat <- ccType(ccGeneList)
	  hasGO <- testCat %in% goType
	  hasKEGG <- testCat %in% keggType
	  if (!(annOpt$godb) && (sum(hasGO) > 0)) {
	    stop('GO requested, but GO.db not loaded. Please load GO.db.', call.=FALSE)
	  } 
	  
	  if (!(annOpt$keggdb) && (sum(hasKEGG) > 0)) {
	    stop('KEGG requested, but KEGG.db not loaded. Please load KEGG.db', call.=FALSE)
	  }
	  
	  anyNames <- sapply(ccGeneList, function(x){
	  	names(x$any.annotation)
	  })
	  anyNames <- paste("ANY.", unique(unlist(anyNames)), sep="")
	  hasANY <- testCat %in% anyNames
	  # get rid of anything else besides GO, KEGG and ANY that is in the requested things
	  testCat <- testCat[(hasGO | hasKEGG | hasANY)]
	  nCat <- length(testCat)
	  nList <- length(ccNames)
	  
	  allEnrich <- vector("list", nCat)
	  
	  names(allEnrich) <- testCat
	  
	  cat("Performing Enrichment Calculations ....\n")
	  
	  for (iCat in 1:nCat){
	    if (hasGO[iCat]){
	      allEnrich[[iCat]] <- .goEnrich(ccGeneList,testCat[iCat])
	    } else if (hasKEGG[iCat]){
	      allEnrich[[iCat]] <- .keggEnrich(ccGeneList)
	    } else if (hasANY[iCat]) {
	    	allEnrich[[iCat]] <- .anyEnrich(ccGeneList, testCat[iCat])
	    }
	  }
	  
	  cat("Done!!\n\n")
	  allEnrich <- new("ccEnrichCollection", allEnrich)
	  return(allEnrich)
})

.goEnrich <- function(ccGeneList, testCat){
  # assume that anything that made it in here has only "GO" tests to do
  testNames <- names(ccGeneList)
  nTest <- length(testNames)
	allAnn <- vector("list", nTest)
	names(allAnn) <- testNames
	
	# set up the default GOHyperGParamsCC object
	testAnn <- new("GOHyperGParamsCC", geneIds=ccGeneList[[1]]$genes, universeGeneIds=ccGeneList[[1]]$universe, 
				annotation=ccGeneList[[1]]$annotation, ontology=testCat, conditional=FALSE, 
				testDirection=testDirection(ccGeneList), fdr=fdr(ccGeneList), pvalueCutoff = pvalueCutoff(ccGeneList))
  
	for (iTest in 1:nTest) {
		
		# just go through and do all of the required calculations. Will use another function to extract and plot the results
    # This way, we, and the user have an object they can requery using different cutoffs if they want
    # basically trying to separate the enrichment calculations from the generation of tables and plots
    

		# keep the user updated on what is going on
		printStr <- paste("GO ", testCat, ': ', testNames[iTest], '\n', sep=" ")
		cat(printStr)
		
    # set the different geneIds, universeIds, and annotation
		geneIds(testAnn) <- ccGeneList[[iTest]]$genes
    testAnn@universeGeneIds <- ccGeneList[[iTest]]$universe
    annotation(testAnn) <- ccGeneList[[iTest]]$annotation
    # calculate the enriched annotations
		tmpAnn <- hyperGTestCC(testAnn)
    allAnn[[iTest]] <- tmpAnn
	}
  
  allAnn <- new("GOccEnrichResult", allAnn, fdr=fdr(ccGeneList), pvalueCutoff=pvalueCutoff(ccGeneList), ontology=testCat, pvalueType=pvalueType(allAnn[[1]]),minCount=0)
	
  return(allAnn)
}


.keggEnrich <- function(ccGeneList){
  testNames <- names(ccGeneList)
  nTest <- length(testNames)
  allAnn <- vector("list", nTest)
  names(allAnn) <- testNames
  
  testAnn <- new("KEGGHyperGParamsCC", geneIds=ccGeneList[[1]]$genes, universeGeneIds=ccGeneList[[1]]$universe,
        annotation=ccGeneList[[1]]$annotation, testDirection=testDirection(ccGeneList), fdr=fdr(ccGeneList), pvalueCutoff=pvalueCutoff(ccGeneList))
        
  for (iTest in 1:nTest){
    printStr <- paste('KEGG: ', testNames[iTest], '\n', sep="")
    cat(printStr)
    
    geneIds(testAnn) <- ccGeneList[[iTest]]$genes
    testAnn@universeGeneIds <- ccGeneList[[iTest]]$universe
    annotation(testAnn) <- ccGeneList[[iTest]]$annotation
    # calculate the enriched annotations
    tmpAnn <- hyperGTestCC(testAnn)
    allAnn[[iTest]] <- tmpAnn
  }
  
 allAnn <- new("KEGGccEnrichResult", allAnn, fdr=fdr(ccGeneList), pvalueCutoff=pvalueCutoff(ccGeneList), pvalueType=pvalueType(allAnn[[1]]),minCount=0)
 return(allAnn)
}

.anyEnrich <- function(ccGeneList, testCat){
	testNames <- names(ccGeneList)
	nTest <- length(testNames)
#	allAnn <- vector("list", nTest)
# 	names(allAnn) <- testNames
	
	anyTestCat <- substr(testCat, 5, nchar(testCat)) # assumes "ANY." is the first, but we want the last bit
	
	allAnn <- lapply(testNames, function(x){
		message(testCat,': ', x)
# 		cat(printStr)
		testAnn <- new("ANYHyperGParamsCC",
									 geneIds=ccGeneList[[x]]$genes,
									 universeGeneIds=ccGeneList[[x]]$universe,
									 annotation=ccGeneList[[x]]$any.annotation[[anyTestCat]][["annotation"]],
									 testDirection=testDirection(ccGeneList),
									 fdr=fdr(ccGeneList),
									 pvalueCutoff=pvalueCutoff(ccGeneList))
		tmpAnn <- hyperGTestCC(testAnn)
		tmpAnn@link <- ccGeneList[[x]]$any.annotation[[anyTestCat]][["link"]]
		tmpAnn@description <- ccGeneList[[x]]$any.annotation[[anyTestCat]][["description"]]
		tmpAnn@ccType <- testCat
		return(tmpAnn)
	})	
	
	names(allAnn) <- testNames
	
	allAnn <- new("ANYccEnrichResult", allAnn, 
								fdr=fdr(ccGeneList), 
								pvalueCutoff=pvalueCutoff(ccGeneList), 
								pvalueType=pvalueType(allAnn[[1]]),
								minCount=0)
	return(allAnn)
}
# this is an example of the Entrez links we want to make:
# http://www.ncbi.nlm.nih.gov/sites/entrez?cmd=search&db=gene&term=230959[uid]&doptcmdl=Full_Report
# Amigo Links: http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=GO:0010165

# creating links to locations on another page:
# source: <a href="pagename.html#linkname">...</a>, destination: <a name="linkname"></a> Note that for a table, the desitination must be inside of a TD tag
