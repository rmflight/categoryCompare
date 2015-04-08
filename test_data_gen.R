# test data generation for v2 enrichment objects

# because I am concentrating on combining objects first, I need some actual
# enrichment data that can be used
#
# this script generates test data that we can subsequently use for generating and
# combining objects

library("affy")
library("hgu95av2.db")
library("genefilter")
library("estrogen")
library("limma")
library("GOstats")
#library("Category")


datadir <- system.file("extdata", package = "estrogen")
pd <- read.AnnotatedDataFrame(file.path(datadir,"estrogen.txt"), 
                              header = TRUE, sep = "", row.names = 1)
currDir <- getwd()
setwd(datadir)
a <- ReadAffy(filenames=rownames(pData(pd)), phenoData = pd, verbose = TRUE)
setwd(currDir)
eData <- rma(a)

e10 <- eData[, eData$time.h == 10]
e10 <- nsFilter(e10, remove.dupEntrez=TRUE, var.filter=FALSE, 
                feature.exclude="^AFFX")$eset

e10$estrogen <- factor(e10$estrogen)
d10 <- model.matrix(~0 + e10$estrogen)
colnames(d10) <- unique(e10$estrogen)
fit10 <- lmFit(e10, d10)
c10 <- makeContrasts(present - absent, levels=d10)
fit10_2 <- contrasts.fit(fit10, c10)
eB10 <- eBayes(fit10_2)
table10 <- topTable(eB10, number=nrow(e10), p.value=1, adjust.method="BH")
table10$Entrez <- unlist(mget(rownames(table10), hgu95av2ENTREZID, ifnotfound=NA))

e48 <- eData[, eData$time.h == 48]
e48 <- nsFilter(e48, remove.dupEntrez=TRUE, var.filter=FALSE, 
                feature.exclude="^AFFX" )$eset

e48$estrogen <- factor(e48$estrogen)
d48 <- model.matrix(~0 + e48$estrogen)
colnames(d48) <- unique(e48$estrogen)
fit48 <- lmFit(e48, d48)
c48 <- makeContrasts(present - absent, levels=d48)
fit48_2 <- contrasts.fit(fit48, c48)
eB48 <- eBayes(fit48_2)
table48 <- topTable(eB48, number=nrow(e48), p.value=1, adjust.method="BH")
table48$Entrez <- unlist(mget(rownames(table48), hgu95av2ENTREZID, ifnotfound=NA))

gUniverse <- unique(union(table10$Entrez, table48$Entrez))

g10 <- unique(table10$Entrez[table10$adj.P.Val < 0.05])
g48 <- unique(table48$Entrez[table48$adj.P.Val < 0.05])

g10_hyper <- new("GOHyperGParams", geneIds = g10,
                 universeGeneIds = gUniverse,
                 annotation = "org.Hs.eg.db",
                 ontology = "BP",
                 pvalueCutoff = 1.0)

g10_enrich <- hyperGTest(g10_hyper)

g48_hyper <- new("GOHyperGParams", geneIds = g48,
                 universeGeneIds = gUniverse,
                 annotation = "org.Hs.eg.db",
                 ontology = "BP",
                 pvalueCutoff = 1.0)
g48_enrich <- hyperGTest(g48_hyper)
g48_counts <- geneCounts(g48_enrich)
g48_pvals <- pvalues(g48_enrich)

bp_annotation <- select(org.Hs.eg.db, keys = gUniverse, keytype = "ENTREZID", columns = "GOALL")
bp_annotation <- bp_annotation[(bp_annotation$ONTOLOGYALL %in% "BP"),]

bp_annotation <- split(bp_annotation$ENTREZID, bp_annotation$GOALL)
bp_annotation <- lapply(bp_annotation, unique)

estrogen_10_enrichment <- list(features = g10,
                               universe = gUniverse,
                               annotation = bp_annotation,
                               enrich = list(pvalues = pvalues(g10_enrich),
                                             odds_ratio = oddsRatios(g10_enrich),
                                             counts = geneCounts(g10_enrich),
                                             stat_names = names(pvalues(g10_enrich))))

estrogen_48_enrichment <- list(features = g48,
                               universe = gUniverse,
                               annotation = bp_annotation,
                               enrich = list(pvalues = pvalues(g48_enrich),
                                             odds_ratio = oddsRatios(g48_enrich),
                                             counts = geneCounts(g48_enrich),
                                             stat_names = names(pvalues(g48_enrich))))

save(estrogen_10_enrichment, estrogen_48_enrichment, file = "test_data.RData")
