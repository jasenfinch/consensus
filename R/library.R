
fileName <- function(mf,db,organism,path = 'construction_library'){
  
  if (length(organism) == 0) {
    organism <- 'none'
  }
  
  file_name <- switch(db,
                      kegg = paste0(path,'/',mf,'_',db,'_',organism,'.rds'),
                      pubchem = paste0(path,'/',mf,'_',db,'.rds')
  )
  
  return(file_name)
}

#' @importFrom readr write_rds

setGeneric('saveConsensus',function(x,path = 'construction_library')
  standardGeneric('saveConsensus'))

setMethod('saveConsensus',signature = 'Consensus',
          function(x,path = 'construction_library'){
            
            if (!dir.exists(path)){
              dir.create(path,recursive = TRUE)
            }
            
            message(str_c('Exporting to ',path))
            
            write_rds(
              x,
              fileName(
                mf(x),
                database(x),
                organism(x),
                path = path
              ))
          })

#' @importFrom readr read_rds

status <- function(items,library_path,threshold = 50){
  
  item_status <- items %>% 
    rowwise() %>% 
    mutate(
      file_name = fileName(MF,database,organism,path = library_path),
      exists = file.exists(file_name),
      status = if (exists) {
        read_rds(file_name) %>% 
          consensus(Adduct,threshold) %>% 
          .$kingdom %>% 
          .[1]
      } else {
        'Unavailable'
      },
      status = status %>%
        replace(
          !. %in% c(
            'No database hits',
            'Unclassified',
            'Unavailable'
          ),
          'Classified'
        ),
      status = factor(status,
                      levels = c(
                        'Classified',
                        'Unclassified',
                        'No database hits',
                        'Unavailable'
                      ))
    )
  
  return(item_status)
}


searchNecessary <- function(search_list,db){
  
  search_list <- search_list %>% 
    dplyr::filter(
      status != 'Classified'
    )
  
  if ('kegg' %in% db){
    search_list <- search_list %>% 
      dplyr::group_by(MF) %>% 
      dplyr::filter(
        'kegg' %in% database
      )
  }
  
  search_list <- search_list %>% 
    filter(
      status == 'Unavailable'
    )
  
  return(search_list)
}
