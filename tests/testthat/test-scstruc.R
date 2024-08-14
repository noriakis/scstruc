library(scran)
library(scstruc)
library(bnlearn)

test_that("Basic function without errors", {
	sce <- mockSCE()
	sce <- logNormCounts(sce)
	included_genes <- sample(row.names(sce), 10)
	expect_error(scstruc(sce, included_genes, changeSymbol=FALSE), NA)
})

