#' @importFrom readr write_rds

setMethod('saveConsensus',signature = 'Consensus',
          function(x,path = '.'){
            
            message(str_c('Exporting to ',path))
            
            if (database(x) == 'kegg') {
              fileName <- str_c(str_c(mf(x),database(x),organism(x),sep = '_'),'.rds')  
            }
            
            if (database(x) == 'pubchem') {
              fileName <- str_c(str_c(mf(x),database(x),sep = '_'),'.rds')
            }
            
            write_rds(x,str_c(path,fileName,sep = '/'))
          })

checkLibrary <- function(path){
  str_c(path,'construction_library',sep = '/') %in% list.dirs(path)
}

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