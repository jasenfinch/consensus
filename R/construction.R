
globalVariables(c('Compound','Consensus (%)','Enzyme','InChI','SMILES','ID','Level','Label','INCHIKEY'))

#' construction
#' @description Build or add to and load a consensus classification library. 
#' @param MFs vector of molecular formulas
#' @param path target file path for classification library 
#' @param db databases to search. Can be either kegg or pubchem.
#' @param organism KEGG organism ID. Ignored if kegg is not specified in db.
#' @param adductRules adduct rules table as returned by mzAnnotation::adducts()
#' @param threshold \% majority threshold for consensus classifications 
#' @export

construction <- function(MFs, path = '.', db = c('kegg','pubchem'), organism = character(), threshold = 50){
  
  mfs <- MFs %>%
    map(~{
      tibble(MF = .,db = db)
    }) %>%
    bind_rows()
  
  libraryPath <- str_c(path,'classification_library',sep = '/')
  
  if (isFALSE(checkLibrary(path))) {
    dir.create(libraryPath)  
    message(str_c('\nCreated structural classifcation library at ',libraryPath))
  } else {
    message(str_c('\nUsing structural classifcation library at ',libraryPath))
    
    classificationLibrary <- loadLibrary(path)
    
    statuses <- classificationLibrary %>%
      map(status) %>%
      bind_rows()
    
    todoMFs <- mfs %>%
      split(.$MF) %>%
      map(~{
        d <- .
        cl <- d %>%
          left_join(classificationLibrary, by = c("MF", "db")) %>%
          # filter(!is.na(kingdom)) %>%
          select(db,kingdom) %>%
          distinct() %>%
          split(.$db) %>%
          map(~{
            kingdoms <- .$kingdom
            if (!identical(kingdoms,'No hits') & !identical(kingdoms,'Unclassified') & !identical(kingdoms,c('No hits','Unclassified'))) {
              return(NULL)
            } else {
              return(.$db[1])
            }
          }) #%>%
        unlist(use.names = FALSE)
      })
    
    
    
  }
  
  
  mfs %>%
    split(1:nrow(.)) %>%
    walk(~{
      dbase <- .$db %>%
        str_split('; ') %>%
        .[[1]]
      
      for (i in dbase){
        consense <- construct(.,db = i,organism = organism,threshold = threshold)
        
        saveConsensus(consense,path = libraryPath) 
        
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
    })
  
  message('\nComplete!')
  
}
