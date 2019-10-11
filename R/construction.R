
globalVariables(c('Compound','Consensus (%)','Enzyme','InChI','SMILES','ID','Level','Label','INCHIKEY'))

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