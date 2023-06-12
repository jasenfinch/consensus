test_that("sankey plots work", {
  
  sankey_plot <- plotSankey(structural_classifications)
  
  expect_s3_class(sankey_plot,'ggplot')
})
