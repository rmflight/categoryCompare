# test creating a new ccCompare method for ANY objects

data("ccData", package="categoryCompare")

library(categoryCompareANY)
library(org.Hs.eg.db)
library(GO.db)

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
allHsPub <- allHsPub[(pubLen <= 500)]

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

enrichList <- ccEnrich(geneList)

ccOpts <- new("ccOptions", listNames=names(geneList))

pvalueType(enrichList) <- "pval" # needed to make sure that comparing BP and ANY BP are the exact same

compareList <- ccCompare(enrichList, ccOpts)

extGO <- categoryCompareANY:::.extractRes(enrichList$BP)
extANYGO <- categoryCompareANY:::.extractRes(enrichList$ANY.GOBP)
extANYGO$allNodes[!(extANYGO$allNodes %in% extGO$allNodes)]
length(extANYGO$allNodes)
length(extGO$allNodes)

extANYGO$sigID[!(extANYGO$sigID %in% extGO$sigID)]
length(extANYGO$sigID)
length(extGO$sigID)

goCatGene <- extGO$allCatGene[(names(extGO$allCatGene) %in% extGO$sigID)]
anyCatGene <- extANYGO$allCatGene[(names(extANYGO$allCatGene) %in% extANYGO$sigID)]

goCatGene <- goCatGene[order(names(goCatGene))]
goCatGene <- lapply(goCatGene, function(x){sort(x)})
anyCatGene <- anyCatGene[order(names(anyCatGene))]
anyCatGene <- lapply(anyCatGene, function(x){sort(x)})

all.equal(goCatGene, anyCatGene) # returns TRUE!!!

useANYGene <- anyCatGene[names(anyCatGene) %in% extANYGO$sigID]
useGOGene <- categoryCompareANY:::getGO2ALLEGS(extGO$sigID, "org.Hs.eg.db")

anyGraph <- categoryCompareANY:::createGraph2(extANYGO$sigID, useANYGene, "GO")
goGraph <- categoryCompareANY:::createGraph2(extGO$sigID, useGOGene, "GO")

anyGraph
goGraph

# this appears to be an odd result. Having **more** of the GO to Entrez annotations actually
# results in fewer links and nodes in the final graph. This was actually rather counter-intuitive
# to what I expected. But the good news is that categorgyCompareANY does appear to work. 

# so why are the graphs different?
tmpBP <- compareList$BP@mainGraph
tmpA <- compareList$ANY.GOBP@mainGraph

goEdges <- edges(tmpBP)
anyEdges <- edges(tmpA)
anyinGO <- anyEdges[(names(anyEdges) %in% names(goEdges))]
anyinGO <- anyinGO[order(names(anyinGO))]
anyinGO <- lapply(anyinGO, sort)

goinANY <- goEdges[(names(goEdges) %in% names(anyEdges))]
goinANY <- goinANY[order(names(goinANY))]
goinANY <- lapply(goinANY, sort)

all.equal(anyinGO, goinANY)

goDat <- edgeData(tmpBP, names(goinANY)[1], goinANY[[1]])
anyDat <- edgeData(tmpA, names(anyinGO)[1], anyinGO[[1]])


goMatch <- which(names(goDat) %in% names(anyDat))
anyMatch <- which(names(anyDat) %in% names(goDat))

goDat[goMatch]
anyDat[anyMatch]

debug(categoryCompareANY:::.ccCompareANY)
categoryCompareANY:::.ccCompareANY(enrichList$ANY.GOBP, ccOpts)
# now why does the ANY.GOBP have more than the BP?

tmpBP <- compareList$BP@mainGraph
tmpA <- compareList$ANY.GOBP@mainGraph

aNodes <- nodes(tmpA)[!(nodes(tmpA) %in% nodes(tmpBP))]