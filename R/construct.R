globalVariables(c('.','kingdom','CID','MF','Adduct','InChIKey','superclass','subclass','level 5','MF',
                  'Charge','CanonicalSMILES','CovalentUnitCount','SMILE','ACCESSION_ID','MolecularFormula','com','Name','INCHI','Score','Feature','Intensity'))

#' @importFrom mzAnnotation adducts
#' @importFrom methods new

construct <- function(MF, db, organism = character(), threshold = 50){
  
  if (!(db %in% c('kegg','pubchem'))) {
    stop('Database not recognised!')
  }
  
  adductRules = adducts()
  
  message(str_c('\n',MF))
  
  consense <- new('Consensus')
  consense@MF <- MF
  consense@database <- db
  consense@adductRules <- adductRules
  consense@threshold <- threshold
  
  if (i == 'kegg') {
    consense@organism <- organism
  }  
  
  if (i == 'pubchem') {
    consense@organism <- character()
    consense@database <- 'pubchem' 
  }
  
  consense <- mfHits(consense)
  consense <- classify(consense)
  consense <- pips(consense)
  consense <- consensus(consense)
  
  return(consense)
}
