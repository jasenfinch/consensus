#' Consensus
#' @description An S4 class to store consensus structral classification results for a molecular formula.
#' @slot MF molecular formula
#' @slot adductRules tibble containing adduct rules as returned by mzAnnotation::adducts()
#' @slot organism organism KEGG ID. NA if database is pubchem.
#' @slot database database, kegg or pubchem
#' @slot threshold % majority for consensus classification
#' @slot hits MetaboliteDatabase containing molecular formula hits from given database
#' @slot classifications structural classifications for database hits
#' @slot PIPs putative ionisation product matches from database hits for given adduct rules
#' @slot consensus consensus structural classifications for adducts
#' @importClassesFrom mzAnnotation MetaboliteDatabase
#' @export

setClass('Consensus',
         slots = list(
           MF = 'character',
           adductRules = 'tbl_df',
           organism = 'character',
           database = 'character',
           threshold = 'numeric',
           hits = 'MetaboliteDatabase',
           classifications = 'tbl_df',
           PIPs = 'tbl_df',
           consensus = 'tbl_df'
         )
)
