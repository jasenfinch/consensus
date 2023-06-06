
test_that('construction works',{
  MFs <- tibble(
    MF = c('C18H29NO14',
           'C17H32O13S',
           'C14H6OS9'),
    Adduct = '[M-H]1-'
  )
  
  ClassyFireDB <- RSQLite::dbConnect(RSQLite::SQLite(), 
                                     paste0(tempdir(),'ClassyFireCache.db'))
  
  structural_classifications <- construction(MFs,
                                             db = c('kegg',
                                                    'pubchem'),
                                             organism = 'bdi',
                                             classyfireR_cache = paste0(tempdir(),'ClassyFireCache.db'))
  
  expect_s3_class(structural_classifications,'tbl_df')
  
  structural_classifications <- construction(MFs,
                                             db = c('kegg',
                                                    'pubchem'))
  
  
  expect_s3_class(structural_classifications,'tbl_df')
  
  unlink(paste0(tempdir(),'ClassyFireCache.db'))
})

test_that('construction works for the Assignment S4 class',{
  p <- assignments::assignmentParameters('FIE-HRMS')
  
  assignments <- assignments::assignMFs(
    assignments::feature_data,p)
  
  assignments@assignments <- assignments@assignments %>% 
    filter(MF == 'C6H9NO7')
  
  structural_classifications <- construction(
    assignments,
    db = c('kegg'))
  
  expect_output(print(structural_classifications),'Consensus')
  expect_s4_class(structural_classifications,'Construction')
  expect_s3_class(classifications(structural_classifications),'tbl_df')
  expect_s3_class(summariseClassifications(structural_classifications),'tbl_df')
})
