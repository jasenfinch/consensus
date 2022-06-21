
context('construction')

test_that('construction workds',{
  MFs <- tibble(MF = c('C4H11N','C10H14O7P2'),
              Adduct = rep('[M-H]1-',2))
  
  res <- construction(MFs)
  res1 <- construction(MFs)
  
  expect_equal(nrow(res),2)
  expect_equal(nrow(res1),2)
})
