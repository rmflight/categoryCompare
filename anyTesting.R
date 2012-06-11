# testing to see if we get the same results on GO BP using categoryCompare and 
# categoryCompareANY

# load up some test data
library(categoryCompare)
library(org.Hs.eg.db)
library(GO.db)
data(ccData)

g48 <- table48$Entrez

hyperParam1 <- new("GOHyperGParamsCC", 
									 geneIds=g48, 
									 universeGeneIds=gUniverse, 
									 categoryName="GO",
									 ontology="BP",
									 annotation="org.Hs.eg.db",
									 fdr=0)

hyperGRes1 <- hyperGTestCC(hyperParam1)


detach(package:categoryCompare)
library(categoryCompareANY)

hyperParam2 <- new("GOHyperGParamsCC", 
									 geneIds=g48, 
									 universeGeneIds=gUniverse, 
									 categoryName="GO",
									 ontology="BP",
									 annotation="org.Hs.eg.db",
									 fdr=0)

hyperGRes2 <- hyperGTestCC(hyperParam2)

# now build a fake annotation based on GO
allHsGO <- as.list(org.Hs.egGO2ALLEGS)

allGOIDs <- names(allHsGO)
goOnt <- Ontology(allGOIDs)
allHsGO <- allHsGO[(goOnt %in% "BP")]

ent2GO <- reverseSplit(allHsGO)
ent2GO <- ent2GO[(names(ent2GO) %in% gUniverse)]

useGO <- reverseSplit(ent2GO)
useGO <- sapply(useGO, unique) # turns out this is critical. *Most* annotations should not be a problem, but probably a good thing to check

hyperParam3 <- new("ANYHyperGParamsCC",
									 geneIds=g48,
									 universeGeneIds=gUniverse,
									 categoryName="ANY",
									 annotation=useGO,
									 fdr=0)

hyperGRes3 <- hyperGTestCC(hyperParam3)

pval3 <- hyperGRes3@pvalues
pval2 <- hyperGRes2@pvalues

pval2 <- pval2[order(names(pval2))]
pval3 <- pval3[order(names(pval3))]

cat2Gene2 <- hyperGRes2@catToGeneId
cat2Gene3 <- hyperGRes3@catToGeneId

cat2Gene2 <- cat2Gene2[order(names(cat2Gene2))]
cat2Gene3 <- cat2Gene3[order(names(cat2Gene3))]

res2Cnt <- sapply(cat2Gene2, length)
res3Cnt <- sapply(cat2Gene3, length)

### now lets think about how we want to set up our annotations to
# enable multiple types of custom annotations, and use them later
library(categoryCompareANY)
library(org.Hs.eg.db)
library(GO.db)
data(ccData)

g48 <- table48$Entrez
g10 <- table10$Entrez

# set up a custom GOBP annotation
allHsGO <- as.list(org.Hs.egGO2ALLEGS)

allGOIDs <- names(allHsGO)
goOnt <- Ontology(allGOIDs)
allHsGO <- allHsGO[(goOnt %in% "BP")]

ent2GO <- reverseSplit(allHsGO)
ent2GO <- ent2GO[(names(ent2GO) %in% gUniverse)]

useGO <- reverseSplit(ent2GO)
useGO <- sapply(useGO, unique)

GOBP <- list(annotation=useGO, description=Term(names(useGO)))

# and how about a PubMed based annotation
allHsPub <- mget(gUniverse, org.Hs.egPMID, ifnotfound=NA) # only those in gUniverse
allHsPub <- reverseSplit(allHsPub)
allHsPub <- sapply(allHsPub, unique)
allHsPub <- sapply(allHsPub, function(x){x[!(is.na(x))]}) # make sure to get rid of NA

# take out anything with more than 10 genes (probably a sequencing project?)
pubLen <- sapply(allHsPub, length)
allHsPub <- allHsPub[(pubLen <= 10)]

pubLink <- paste("http://www.ncbi.nlm.nih.gov/pubmed/", names(allHsPub), sep="")
names(pubLink) <- names(allHsPub)
PUB <- list(annotation=allHsPub, link=pubLink)

geneList <- list(g48=list(genes=g48,
													universe=gUniverse,
													annotation="org.Hs.eg.db",
													any.annotation=list(GOBP=GOBP, PUB=PUB)),
								 g10=list(genes=g10,
								 				 universe=gUniverse,
								 				 annotation="org.Hs.eg.db",
								 				 any.annotation=list(GOBP=GOBP, PUB=PUB)))

geneList <- new("ccGeneList", geneList, ccType=c("BP", "ANY.GOBP", "ANY.PUB"))
fdr(geneList) <- 0
enrichRes <- ccEnrich(geneList)
pvalueCutoff(enrichRes) <- 0.001
ccOptions <- new("ccOptions", listNames=names(geneList))
compareRes <- ccCompare(enrichRes, ccOptions)

compareRes$ANY.PUB <- breakEdges(compareRes$ANY.PUB, 0.8)
cw.PUB <- ccOutCyt(compareRes$ANY.PUB, ccOptions)

# what is going on
x <- geneList[[1]]

testAnn <- new("ANYHyperGParamsCC",
							 geneIds=x$genes,
							 universeGeneIds=x$universe,
							 annotation=x$any.annotation[[anyTestCat]][["annotation"]],
							 testDirection=testDirection(geneList),
							 fdr=fdr(geneList),
							 pvalueCutoff=pvalueCutoff(geneList))
tmpAnn <- hyperGTestCC(testAnn)

testGO <- new("GOHyperGParamsCC",
							geneIds=x$genes,
							universeGeneIds=x$universe,
							annotation=x$annotation,
							ontology="BP",
							fdr=fdr(geneList),
							pvalueCutoff=pvalueCutoff(geneList))
tmpGO <- hyperGTestCC(testGO)


# need to get ALL of the category to geneIds into a single list
tmpRes <- enrichRes[[2]]
tmpCatGeneNames <- unique(unlist(sapply(tmpRes, function(x){names(x@catToGeneId)})))
allCatGene <- vector('list', length(tmpCatGeneNames))
names(allCatGene) <- tmpCatGeneNames

nList <- length(tmpRes)

addGene <- function(catName){
	unique(c(allCatGene[[catName]], tmpCatGene[[catName]]))
}

for (iList in 1:nList){
	tmpCatGene <- tmpRes[[iList]]@catToGeneId
	allCatGene <- sapply(tmpCatGeneNames, addGene)
}