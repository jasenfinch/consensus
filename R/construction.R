
globalVariables(c('Compound','Consensus (%)','Enzyme','InChI','SMILES','ID','Level','Label','INCHIKEY'))

#' construction
#' @description Build or add to a consensus classification library. 
#' @param MFs vector of molecular formulas
#' @param path target file path for classification library 
#' @param db databases to search. Can be either kegg or pubchem.
#' @param organism KEGG organism ID. Ignored if kegg is not specified in db.
#' @param adductRules adduct rules table as returned by mzAnnotation::adducts()
#' @param threshold \% majority threshold for consensus classifications 
#' @export

construction <- function(MFs, path = '.', db = c('kegg','pubchem'), organism = character(), adductRules = adducts(), threshold = 50){
  
  libraryPath <- str_c(path,'construction_library',sep = '/')
  
  if (isFALSE(checkLibrary(path))) {
    dir.create(libraryPath)  
  }
  
  MFs %>%
    map(~{
      construct(.,db = db,organism = organism,adductRules = adductRules,threshold = threshold) %>%
        saveConsensus(path = libraryPath)
    })
  
  message('\nComplete!')
  
}