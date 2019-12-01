#test_DGOEnrichmentAnalysis.R



context("DGOEnrichmentAnalysis")



# ==== LOAD DATA =================================================

#

data(geneList)  # load data from DOSE package


#

# ==== END SETUP AND PREPARE ===================================================


gene <- names(geneList)[abs(geneList) > 2]
DGOResult <- enrichDGO(gene, universe = names(geneList))

DOanalysis <- DGOResult[["DO"]]
GOanalysis <- DGOResult[["GO"]]

test_that("corrupt input generates errors",  {
  
  expect_error(enrichDGO(c(), universe = names(geneList)),
               "The input should be a vector of gene IDs with a fold change")
  
  
})



test_that("a sample input prouces the expected output",  {
  
  expect_equal(class(DOanalysis), class(GOanalysis))
  
})





# ==== BEGIN TEARDOWN AND RESTORE ==============================================

# Remove any variables that the test has created 

#

rm(DOanalysis)

rm(GOanalysis)




# ==== END  TEARDOWN AND RESTORE ===============================================



# [END]