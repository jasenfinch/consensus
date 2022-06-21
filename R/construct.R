globalVariables(c('.','kingdom','CID','MF','Adduct','InChIKey','superclass','subclass','level 5','MF',
                  'Charge','CanonicalSMILES','CovalentUnitCount','SMILE','id','MolecularFormula','com','Name','INCHI','Score','Feature','Intensity'))

#' @importFrom mzAnnotation adduct_rules
#' @importFrom methods new

construct <- function(MF,db,organism = character(),threshold = 50,adduct_rules_table = adduct_rules()){
  
  if (!(db %in% c('kegg','pubchem'))) {
    stop('Database not recognised!')
  }
  
  message(str_c('\n',MF))
  
  consense <- new('Consensus')
  consense@MF <- MF
  consense@database <- db
  consense@adductRules <- adduct_rules
  consense@threshold <- threshold
  
  if (db == 'kegg') {
    consense@organism <- organism
  }  
  
  if (db == 'pubchem') {
    consense@organism <- character()
    consense@database <- 'pubchem' 
  }
  
  consense %>%
    mfHits() %>%
    classify() %>%
    pips() %>%
    consensus()
  
}
