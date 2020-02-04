
context('construction')

test_that('construction workds',{
  a <- tibble(MF = c('C4H11N','C10H14O7P2'),Adduct = rep('[M-H]1-',2))  
  res <- construction(a,path = tempdir())
  res1 <- construction(a,path = tempdir())
  
  l <- loadLibrary(a,tempdir())
  
  p <- capture.output(print(l[[1]]))
  
  expect_equal(nrow(res),2)
  expect_equal(nrow(res1),2)
  expect_length(p,13)
})
