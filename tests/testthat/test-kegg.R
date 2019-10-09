
context('KEGG searches')

test_that('KEGG compound info returned correctly',{
  compounds <- keggCompoundInfo(c('C00089','C00149'))
  
  expect_true(identical(class(compounds),c('tbl_df','tbl','data.frame')))
  expect_equal(nrow(compounds),2)
  expect_equal(ncol(compounds),5)
})

test_that('KEGG consensus returned',{
  d <- tibble(MF = c('C4H6O5'),Adduct = '[M-H]1-')
  
  consensi <- keggConsensus(d,'hsa')
  
  expect_true(class(consensi) == 'Consensus')
})