#' @importFrom readr write_rds

setGeneric('saveConsensus',function(x,path = 'structural_classification_library'){
  standardGeneric('saveConsensus')
})

setMethod('saveConsensus',signature = 'Consensus',
          function(x,path = 'structural_classification_library'){
            
            message(str_c('Exporting to ',path))
            
            if (database(x) == 'kegg') {
              fileName <- str_c(str_c(mf(x),database(x),organism(x),sep = '_'),'.rds')  
            }
            
            if (database(x) == 'pubchem') {
              fileName <- str_c(str_c(mf(x),database(x),sep = '_'),'.rds')
            }
            
            write_rds(x,str_c(path,fileName,sep = '/'))
          })

#' @importFrom readr read_rds

loadLibrary <- function(MFs, path = 'structural_classification_library'){
  
  if (!dir.exists(path)) {
    dir.create(path,recursive = TRUE)  
    message(str_c('Created structural classifcation library at ',path))
  } else {
    message(str_c('Using structural classifcation library at ',path))
  }
  
  libraryContents <- list.files(path,full.names = TRUE,pattern = '.rds') %>%
    tibble(MF = basename(.) %>%
             str_split_fixed('_',2) %>%
             .[,1],
           path = .) %>%
    filter(MF %in% MFs$MF)
  
  libraryContents %>%
    .$path %>%
    map(read_rds)
}

setGeneric('status',function(x){
  standardGeneric('status')
})

setMethod('status',signature = 'Consensus',
          function(x){
            
            if (length(organism(x)) == 0) {
              org <- 'none'
            } else {
              org <- organism(x)
            }
            
            st <- x %>%
              consensusClassifications() %>%
              .$kingdom %>%
              unique() %>%
              sort()
            
            if (length(st) == 1 & !identical(st,'No hits') & !identical(st,'Unclassified')) {
              st <- 'Classified'
            }
            
            if (length(st) > 1) {
              
              if ('Unclassified' %in% st) {
                if (identical(st,c('No hits','Unclassified'))) {
                  st <- 'Unclassified'
                } else {
                  st <- 'Classified'
                }  
              } else {
                st <- 'Classified'
              }
            }
            
            
            tibble(MF = mf(x),organism = org,database = database(x),status = st)
          }
)