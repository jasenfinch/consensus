
context('KEGG searches')

test_that('KEGG compound info returned correctly',{
  compounds <- keggCompoundInfo(c('C00089','C00149'))
  
  expect_true(identical(class(compounds),c('tbl_df','tbl','data.frame')))
  expect_equal(nrow(compounds),2)
  expect_equal(ncol(compounds),5)
})