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



test_that("corrupt input generates errors", {
  
  expect_error(DGObarplot(c()),
               "The input should be a list of 2 enricResult objects.")
  
  r <- c(1,2)
  expect_error(DGObarplot(r),
               "DO enrichment results and GO enrichment results")
  
  names(r) <- c("DO", "GO")
  expect_error(DGObarplot(r), "2 of enrichResult objects")
  
})

test_that("valid number of gene ontology groups and showCategory value", {
  
  # check warning message properly produces when
  # showCategory is greater than the number of gene ontology groups
  # in DGO enrichment analysis
  expect_warning(DGObarplot(DGOResult, showCategory = 8),
                 "Input showCategory")
  
  # make an empty copy of DGOResult to test that DGObarplot throws an error
  DGObarEmpty <- DGOResult
  DGObarEmpty[["DO"]]@result <- data.frame()
  DGObarEmpty[["GO"]]@result <- data.frame()
  
  expect_error(DGObarplot(DGObarEmpty, showCategory = 8),
               "0 ontology groups")

})

test_that("a sample input produces the expected output",  {
  
  expect_equal(class(DGObargraph), c("gg", "ggplot"))
  
})


# [END]