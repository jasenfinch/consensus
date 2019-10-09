#' Consensus
#' @description An S4 class to store consensus structral classification results for a molecular formula.
#' @slot MF molecular formula
#' @slot adductRules tibble containing adduct rules as returned by mzAnnotation::adducts()
#' @slot organism organism KEGG ID. NA if database is pubchem.
#' @slot database database, kegg or pubchem
#' @slot threshold % majority for consensus classification
#' @slot hits molecular formula hits from given database
#' @slot PIPs putative ionisation product matches from database hits for given adduct rules
#' @slot classifications structural classifications for database hits 
#' @slot consensus consensus structural classifications for adducts
#' @export

setClass('Consensus',
         slots = list(
           MF = 'character',
           adductRules = 'tbl_df',
           organism = 'character',
           database = 'character',
           threshold = 'numeric',
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