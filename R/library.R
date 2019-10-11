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