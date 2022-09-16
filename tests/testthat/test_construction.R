
test_that('construction works',{
  MFs <- tibble(
    MF = c('C18H29NO14','C17H32O13S'),
              Adduct = '[M-H]1-'
    )

  structural_classifications <- construction(MFs,
                                             db = c('kegg',
                                                    'pubchem'))

  expect_s3_class(structural_classifications,'tbl_df')
})
