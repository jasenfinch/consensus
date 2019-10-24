
context('construct')

test_that('construct works',{
  a <- construct('C4H11N',db = 'kegg',organism = 'bdi')
  
  expect_true(class(a) == 'Consensus')
})