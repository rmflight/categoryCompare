context("ccListsMethodsTests")

require(org.Hs.eg.db)
require(GO.db)
require(KEGG.db)
dataEnv <- new.env()
data(ccData, package="categoryCompare", envir=dataEnv)
g10 <- unique(dataEnv$table10$Entrez)
g48 <- unique(dataEnv$table48$Entrez)

test_that("can create new ccGeneList objects", {
	list10 <- list(genes=g10, universe=gUniverse, annotation='org.Hs.eg.db')
	list48 <- list(genes=g48, universe=gUniverse, annotation='org.Hs.eg.db')
	
	geneLists <- list(T10=list10, T48=list48)
	geneLists <- new("ccGeneList", geneLists, ccType="CC")
	
})


rm(dataEnv)