
context('construction')

test_that('construction workds',{
  a <- tibble(MF = c('C4H11N','C10H14O7P2'),Adduct = rep('[M-H]1-',2))  
  res <- construction(a,path = tempdir())
  
  expect_equal(nrow(res),2)
})
