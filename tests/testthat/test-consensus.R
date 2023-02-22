test_that("consensus print method works", {
  expect_output(new('Consensus') %>% 
                  print(),
                'Consensus')
  
  expect_output(construct('C18H29NO14',
                          db = 'pubchem') %>% 
                  print(),
                'Consensus')
})
