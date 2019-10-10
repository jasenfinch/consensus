globalVariables(c('.','kingdom','CID','MF','Adduct','InChIKey','superclass','subclass','level 5','MF',
                  'Charge','CanonicalSMILES','CovalentUnitCount','SMILE','ACCESSION_ID','MolecularFormula','com','Name','INCHI','Score','Feature','Intensity'))

#' construct
#' @description Construct consensus structural classifications for a given molecular formula.
#' @param MF molecular formula
#' @param db databases to search. Can be either kegg or pubchem.
#' @param organism KEGG organism ID. Ignored if kegg is not specified in db.
#' @param adductRules adduct rules table as returned by mzAnnotation::adducts()
#' @param threshold % majority threshold for consensus classifications 
#' @importFrom mzAnnotation adducts
#' @export

construct <- function(MF, db = c('kegg','pubchem'), organism = character(), adductRules = adducts(), threshold = 50){
  consense <- new('Consensus')
  consense@MF <- MF
  consense@adductRules <- adductRules
  consense@organism <- organism 
  consense@threshold <- threshold

  if ('kegg' %in% db) {
    consense@database <- 'kegg' 
  }
  # consense@database <- 'pubchem'
  
  
  consense <- mfHits(consense)
  consense <- classify(consense)
  consense <- pips(consense)
  consense <- consensus(consense)
  
  return(consense)
}