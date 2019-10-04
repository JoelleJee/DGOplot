#test_DGOdoubleBarGraph.R



context("DGOdoubleBarGraph")



# ==== LOAD DATA =================================================

#


data(geneList)

gene <- names(geneList)[abs(geneList) > 2]

DGOResult <- enrichDGO(gene, universe = names(geneList))

DGObargraph <- DGObarplot(DGOResult)



#

# ==== END SETUP AND PREPARE ===================================================



test_that("corrupt input generates errors",  {
  
  expect_error(DGObarplot(c()),
               "The input should be a list of 2 enricResult objects.")
  
})



test_that("a sample input prouces the expected output",  {
  
  expect_equal(class(DGObargraph), c("gg", "ggplot"))
  
})





# ==== BEGIN TEARDOWN AND RESTORE ==============================================

# Remove any variables that the test has created 

#



# ==== END  TEARDOWN AND RESTORE ===============================================



# [END]