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
  str_c(path,'classification_library',sep = '/') %in% list.dirs(path)
}

#' @importFrom readr read_rds

loadLibrary <- function(path = '.'){
  libraryPath <- str_c(path,'classification_library',sep = '/')
  
  message(str_c('Loading structural classification library at ',libraryPath))
  
  libraryContents <- list.files(libraryPath,full.names = TRUE)
  
  classificationLibrary <- libraryContents %>%
    map(~{
      consense <- read_rds(.) 
      
      org <- organism(consense)
      
      if (length(org) == 0) {
        org <- NA
      }
      
      consense %>%
        consensusClassifications() %>%
        mutate(db = database(consense),
                organism = org,
                MF = mf(consense)) %>%
        select(db:MF,everything())
      }) %>%
    bind_rows()
  
  return(classificationLibrary)
}

setMethod('status',signature = 'Consensus',
          function(x){
            st <- x %>%
              consensusClassifications() %>%
              .$kingdom %>%
              unique() %>%
              sort()
            
            if ((length(st) == 1 & !identical(st,'No hits') & !identical(st,'Unclassified')) | !identical(st,c('No hits','Unclassified'))) {
              st <- 'Classified'
            }
            
            if (identical(st,c('No hits','Unclassified'))) {
              st <- 'Unclassified'
            }
            
            tibble(MF = mf(x),organism = organism(x),database = database(x),status = st)
          }
)