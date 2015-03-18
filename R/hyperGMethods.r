setMethod("hyperGTestCC",
  signature(p="HyperGParamsCC"),
  function(p) hyperGTestCat(p))

# this function is essentially a modified version of .hyperGTestInternal from the package Category, with the difference that it will do replicate series of calculations using random subsets from the universe list to provide an estimate of the false discovery rate (fdr)

hyperGTestCat <- function(p, className="HyperGResultCC") {
    #browser(text="", condition=NULL, expr=TRUE)
	className <- paste(categoryName(p), className, sep="")
  p <- makeValidParams(p)
  origGeneIds <- geneIds(p)
  p@universeGeneIds <- universeBuilder(p)
  selected <- intersect(geneIds(p), universeGeneIds(p))
  # I know this is cheating, but I didn't want to debug the categoryToEntrezBuilder-methods code and chance breaking something
	geneIds(p) <- p@universeGeneIds 
  cat2Entrez <- categoryToEntrezBuilder(p)
	# assuming that everything has been returned, not just what we are interested in by the source
	useCat2Entrez <- whichCat(p,cat2Entrez,selected)
	
	nFDR <- fdr(p)
  stats <- .doHyperGTest(p, useCat2Entrez, list(),
                           selected)

	ord <- order(stats$p)
	pType <- "pval"
	# if the user has asked for FDR runs to be done
	if (nFDR > 0){
		pType <- "fdr"
		nGene <- length(selected)
		tmpPval <- stats$p[ord]
		loCnt <- matrix(0, nrow=length(tmpPval), ncol=nFDR)
		
		# do the FDR calculations
		for (i in 1:nFDR) {
			# take a random sample from the gene universe and do the exact same calculation we did above
			tmpSelected <- sample(universeGeneIds(p), nGene)
			tmpCat2Entrez <- whichCat(p,cat2Entrez,tmpSelected)
			tmpStats <- .doHyperGTest(p, tmpCat2Entrez, list(), tmpSelected)
			tmpP <- tmpStats$p
			tmpCnt <- sapply(tmpPval, function(x) {
				sum(tmpP <= x)
			})
			loCnt[,i] <- tmpCnt
		}
		
		mnCnt <- apply(loCnt, 1, mean)
		orgCnt <- sapply(tmpPval, function(x) {
			sum(tmpPval <= x)
		})
		# FDR is the average number of times that we found lower pvalues divided by the number of items with an actual lower pvalue
		fdrVal <- as.numeric(mnCnt / orgCnt)
	} else { fdrVal <- numeric(length=length(stats$p)) + 1 }
		
		# put the identifiers on, and set anything greater than 1 to 1
		names(fdrVal) <- names(stats$p[ord])
		fdrVal[fdrVal > 1] <- 1
	# now spit out our results
		new(className,
		pvalues=stats$p[ord],
		oddsRatios=stats$odds[ord],
		expectedCounts=stats$expected[ord],
		catToGeneId=useCat2Entrez[ord],
		fdr=nFDR,
		fdrvalues=fdrVal,
		annotation=annotation(p),
		geneIds=selected,
		testName=categoryName(p),
		pvalueCutoff=pvalueCutoff(p),
		testDirection=testDirection(p),
		organism=organism(p),
		pvalueType=pType,
		minCount=0)
}

whichCat <- function(p,cat2Entrez,geneIds) {
	keep.all <- switch(testDirection(p),
					over=FALSE,
					under=TRUE,
					stop("Bad testDirection slot"))
	
	if (!keep.all){
		# this seems to be slower, but avoids any chance of worrying about the length of list IDs
		numFound <- sapply(cat2Entrez,function(x) sum(geneIds %in% x))
		keepCat <- numFound > 0
		cat2Entrez[keepCat]
		
	} else { cat2Entrez }

}

