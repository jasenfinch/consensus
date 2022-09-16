
context('PubChem searches')

test_that('Pubchem MF matching works',{
  mf <- pubchemMatch('C18H29NO14')
  
  expect_s3_class(mf,'tbl_df')
})
  