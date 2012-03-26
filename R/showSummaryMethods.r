setMethod("summary", signature(object="HyperGResultCC"),
  		function(object, pvalue=pvalueCutoff(object), pType=object@pvalueType, minCount=object@minCount)
			{
				wanted <- isResult(object, pvalue, pType, minCount)
				pvals <- pvalues(object)
				fvals <- fdrvalues(object)
				if (!any(wanted)) {
					warning("No results met the specified criteria.  ",
							"Returning 0-row data.frame", call.=FALSE)
					catIds <- character(0)
					pvals <- odds <- ecounts <- fvals <- numeric(0)
					counts <- ucounts <- integer(0)
				} else {
					pvals <- pvals[wanted]
					fvals <- fvals[wanted]
					ucounts <- universeCounts(object)[wanted]
					catIds <- names(pvals)
					odds <- oddsRatios(object)[wanted]
					ecounts <- expectedCounts(object)[wanted]
					counts <- geneCounts(object)[wanted]
				}
				df <- data.frame(ID=catIds, Pvalue=pvals, FDR=fvals, 
				OddsRatio=odds, ExpCount=ecounts, Count=counts, 
				Size=ucounts, stringsAsFactors=FALSE, row.names=NULL)
				df
			})
		
setMethod("show", signature(object="HyperGResultCC"),
        function(object) {
          cat(description(object), "\n")
          isSig <- (pCC(object,pvalueType(object)) <= object@pvalueCutoff[1])
          isCount <- (geneCounts(object) >= object@minCount)
          nPass <- sum(isSig & isCount)
          cat(length(pvalues(object)), testName(object), "ids tested ")
          cat("(", nPass, " have p <= ", object@pvalueCutoff[1],
          " & count >= ", minCount(object), ")\n", sep="")
          cat("Selected gene set size:", geneMappedCount(object), "\n")
          cat("    Gene universe size:", universeMappedCount(object), "\n")
          cat("    Annotation package:", annotation(object), "\n")
})

setMethod("show", signature(object="ccOptions"),
      function(object) {
        cat("List Names: ", paste(listNames(object), collapse="; "), "\n")
        cat("Comparisons: ", paste(compareNames(object), collapse="; "), "\n")
        cat("Colors: ", paste(compareColors(object), collapse="; "), "\n")
        cat("Output Types: ", paste(outType(object), collapse="; "), "\n")
      })
      
setMethod("show", signature(object="ccGeneList"),
      function(object) {
        allNames <- names(object)
        nLists <- length(allNames)
        for (iList in 1:nLists){
          cat("List: ", allNames[iList], "\n")
          cat("Size of gene list: ", as.character(length(object[[iList]]$genes)), "\n")
          cat("Size of gene universe: ", as.character(length(object[[iList]]$universe)), "\n")
          cat("Annotation: ", object[[iList]]$annotation, "\n\n")
        }
        cat("Types of annotations to examine: ", paste(object@ccType, collapse="; "), "\n")
        cat("Number of FDR runs to perform: ", as.character(object@fdr), "\n")
        cat("pValue Cutoff to decide significantly enriched annotations: ", as.character(object@pvalueCutoff), "\n")
        cat("Testdirection: ", object@testDirection, "represented\n")
        
      })
      
      
setMethod("show", signature(object="ccEnrichResult"),
      function(object) {
        subType <- object@ontology
        cat("    Annotation category: ", object@categoryName, " ", subType, "\n")
        cat("               FDR runs: ", object@fdr, "\n")
        cat("Default p-values to use: ", object@pvalueType, "\n")
        cat("                pCutoff: ", object@pvalueCutoff, "\n\n")
        subNames <- names(object)
        nSub <- length(subNames)
        for (iSub in 1:nSub){
          cat("List: ", subNames[iSub], "\n")
          show(object[[iSub]])
          cat("\n")
        }
      })
        
setMethod("show", signature(object="ccEnrichCollection"),
      function(object) {
        nObj <- length(object)
        for (iObj in 1:nObj){
          show(object[[iObj]])
          cat("\n")
        }
      })
	 
  
setMethod("mainGraph", signature(object="ccCompareResult"),
					function(object){
						object@mainGraph
					})

setMethod("mainTable", signature(object="ccCompareResult"),
					function(object){
						object@mainTable
					})

setMethod("allAnnotation", signature(object="ccCompareResult"),
					function(object){
						object@allAnnotation
					})

setMethod("show", signature(object="ccCompareCollection"),
      function(object) {
      	cat("ccCompare results for:\n\n")
        nObj <- length(object)
        for (iObj in 1:nObj){
          show(object[[iObj]])
          cat("\n")
        }
      })

setMethod("show", signature(object="ccCompareResult"),
					function(object){
						subType <- object@ontology
						tmpTable <- mainTable(object)
						cat("Annotation category: ", object@categoryName, " ", subType, "\n")
						cat("Main graph: ")
						show(mainGraph(object))
					}
						)