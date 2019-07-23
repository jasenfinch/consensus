
setClass('Consensus',
         slots = list(
           hits = 'tbl_df',
           PIPs = 'tbl_df',
           classifications = 'tbl_df',
           consensus = 'tbl_df'
         )
)

#' Consensuses
#' @description An S4 class to store consensus structural classifications.
#' @slot consensuses list containing consensus results for each MF
#' @slot results tibble containing consensus structural classifications for each feature
#' @export

setClass('Consensuses',
         slots = list(
           consensuses = 'list',
           results = 'tbl_df'
         )
)