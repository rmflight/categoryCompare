context("hyperGeometricTests")

require(org.Hs.eg.db)
require(GO.db)
data(ccData)
test_that("can calculate basic hypergeometric enrichment", {
	g10 <- unique(table10$Entrez)
	testGO <- new("GOHyperGParamsCC", geneIds=g10, universeGeneIds=gUniverse, 
								annotation="org.Hs.eg.db", ontology="CC", conditional=FALSE, 
								testDirection="over",fdr=0, pvalueCutoff = 0.05)
	resultGO <- hyperGTestCC(testGO)
	resultOld <- enrichLists$CC$T10
	expect_that(resultGO, is_identical_to(resultOld))
})