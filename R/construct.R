#' Consensus structural classifications for a molecular formula
#' @description Generate a consensus structural classification for a molecular formula.
#' @param MF molecular formula to search
#' @param db databases to search. Can be either `kegg` or `pubchem`.
#' @param organism KEGG organism ID to search. Only relevant when argument `db` includes `kegg`
#' @param threshold percentage majority threshold for consensus classification
#' @param adduct_rules_table dataframe containing adduct formation rules. The defaults is `mzAnnotation::adduct_rules()`.
#' @return An S4 object of class `Consensus`.
#' @examples 
#' construct('C4H6O5')
#' @importFrom mzAnnotation adduct_rules
#' @importFrom methods new
#' @export

construct <- function(MF,
                      db = c('kegg','pubchem'),
                      organism = character(),
                      threshold = 50,
                      adduct_rules_table = adduct_rules()){
  
  db <- match.arg(
    db,
    c('kegg','pubchem'))
  
  if (db == 'pubchem') {
    organism <- character()
  }
  
  message(str_c('\n',MF))
  
  x <- new('Consensus',
           MF = MF,
           database = db,
           organism = organism,
           adduct_rules = adduct_rules_table,
           threshold = threshold)
  
  x %>%
    mfHits() %>%
    classify() %>%
    pips() %>%
    consensus()
  
}
