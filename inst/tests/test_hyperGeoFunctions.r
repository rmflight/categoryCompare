context("hyperGeometricTests")

require(org.Hs.eg.db)
require(GO.db)
dataEnv <- new.env()
data(ccData, package="categoryCompare", envir=dataEnv)
g10 <- unique(dataEnv$table10$Entrez)

test_that("can create new HyperGParamsCC objects", {
	testGO <- new("GOHyperGParamsCC", geneIds=g10, universeGeneIds=dataEnv$gUniverse, 
								annotation="org.Hs.eg.db", ontology="CC", conditional=FALSE, 
								testDirection="over",fdr=0, pvalueCutoff = 0.05)
	expect_that(testGO, is_equivalent_to(dataEnv$hyperGParamT))
	testKEGG <- new("KEGGHyperGParamsCC", geneIds=g10, universeGeneIds=dataEnv$gUniverse,
									annotation="org.Hs.eg.db", testDirection="over",
									fdr=0, pvalueCutoff=0.05)
	expect_that(testKEGG, is_a("KEGGHyperGParamsCC"))
	
	expect_that(new("GOHyperGParams", geneIds=g10, universeGeneIds=dataEnv$gUniverse, 
									annotation="org.Hs.eg.db", ontology="CC", conditional=FALSE, 
									testDirection="over",fdr=0), throws_error())
	
})

test_that("can calculate basic hypergeometric enrichment", {
	resultGO <- hyperGTestCC(testGO)
	resultOld <- dataEnv$enrichLists$CC$T10
	expect_that(resultGO, is_identical_to(resultOld))
	expect_that(resultGO@fdrvalues[1], is_equivalent_to(1))
})

test_that("can perform fdr runs", {
	fdr(testGO) <- 50
	resultGO <- hyperGTestCC(testGO)
	expect_that(resultGO@fdrvalues[1], is_equivalent_to(0))
})

rm(dataEnv)