globalVariables(c('.','kingdom','CID','MF','Adduct','InChIKey','superclass','subclass','level 5','MF',
                  'Charge','CanonicalSMILES','CovalentUnitCount','SMILE','ACCESSION_ID','MolecularFormula','com','Name','INCHI','Score','Feature','Intensity'))

#' construct
#' @description Construct consensus structural classifications for a given molecular formula.
#' @param MF molecular formula
#' @param db databases to search. Can be either kegg or pubchem.
#' @param organism KEGG organism ID. Ignored if kegg is not specified in db.
#' @param adductRules adduct rules table as returned by mzAnnotation::adducts()
#' @param threshold \% majority threshold for consensus classifications 
#' @importFrom mzAnnotation adducts
#' @export

construct <- function(MF, db = c('kegg','pubchem'), organism = character(), adductRules = adducts(), threshold = 50){
  
  if (F %in% (db %in% c('kegg','pubchem'))) {
    stop('Database not recognised!')
  }
  
  message(str_c('\n',MF))
  
  consense <- new('Consensus')
  consense@MF <- MF
  consense@adductRules <- adductRules
  consense@threshold <- threshold

  db <- sort(db)
  
  for (i in db) {
    if (i == 'kegg') {
      consense@organism <- organism
      consense@database <- 'kegg' 
    }  
    
    if (i == 'pubchem') {
      consense@organism <- character()
      consense@database <- 'pubchem' 
    }
    
    consense <- mfHits(consense)
    consense <- classify(consense)
    consense <- pips(consense)
    consense <- consensus(consense)
    
    if (database(consense) == 'kegg') {
      kingdoms <- consense %>%
        consensusClassifications() %>%
        .$kingdom %>%
        unique() %>%
        sort()
      
      if (!identical(kingdoms,'No hits') & !identical(kingdoms,'Unclassified') & !identical(kingdoms,c('No hits','Unclassified'))) {
        break()
      }
    }
  }
  
  return(consense)
}
