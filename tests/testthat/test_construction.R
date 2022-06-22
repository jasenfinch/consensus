
# test_that('construction works',{
#   MFs <- tibble(MF = c('C4H11N','C10H14O7P2'),
#               Adduct = rep('[M-H]1-',2))
# 
#   structural_classifications <- construction(MFs,library_path = '.')
#   structural_classifications_using_cache <- construction(MFs,library_path = '.')
# 
#   expect_identical(structural_classifications,
#                    structural_classifications_using_cache)
# })
