
context('PubChem searches')

test_that('Pubchem MF matching works',{
  mf <- pubchemMatch('C12H22O11')
  
  expect_identical(class(mf),c('tbl_df','tbl','data.frame'))
})