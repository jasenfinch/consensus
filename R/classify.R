
#' @importFrom magrittr set_names %>%
#' @importFrom dplyr distinct select bind_rows filter left_join
#' @importFrom classyfireR get_classification classification open_cache
#' @importFrom purrr map_dbl map
#' @importFrom tidyr spread
#' @importFrom tidyselect last_col
#' @importFrom progress progress_bar
#' @importFrom stringr str_c
#' @importFrom tibble is_tibble
#' @importFrom RSQLite dbDisconnect

setGeneric('classify',function(x,classyfireR_cache = NULL){
  standardGeneric('classify')
})

setMethod('classify',signature = 'Consensus',
          function(x,classyfireR_cache = NULL){
            
            if (!is.null(classyfireR_cache))
              classyfireR_cache <- open_cache(dbname = classyfireR_cache)
            
            inchikey <- x %>%
              entries() %>%
              .$INCHIKEY %>% 
              unique()
            
            if (length(inchikey) > 0) {
              message(str_c('Retrieving classifications for ',length(inchikey),' InChIKeys...'))
              
              pb <- progress_bar$new(
                format = "[:bar] :percent eta: :eta",
                total = length(inchikey), clear = FALSE)
              pb$tick(0)
              
              classi <- inchikey %>%
                map(~{
                  out <- capture.output(cl <- get_classification(.x,conn = classyfireR_cache),
                                        type = 'message')
                  
                  if (is.null(cl)) 
                    cl <- tibble(Level = 'kingdom','Classification' = 'Unclassified',CHEMONT = NA) 
                  else {
                    cl <- cl %>%
                      classification()
                    
                    if (length(cl) == 0) {
                      cl <- tibble(Level = 'kingdom','Classification' = 'Unclassified',CHEMONT = NA)
                    }
                  }
                    
                  
                  if (length(out) > 0)
                    if (!grepl('cached',out)) 
                      Sys.sleep(5)
                  
                  pb$tick()
                  return(cl)
                }) %>%
                set_names(inchikey)
              
              if (!is.null(classyfireR_cache))
                dbDisconnect(classyfireR_cache)
              
              classes <- c('kingdom','superclass','class','subclass')
              
              classi <- classi %>%
                .[map_dbl(.,nrow) > 0] %>%
                map(~{
                  .x %>%
                    select(-CHEMONT) %>%
                    spread(.,Level,Classification) 
                }) %>%
                bind_rows(.id = 'INCHIKEY') %>%
                distinct() %>%
                filter(!is.na(kingdom))
              
              if (nrow(classi) > 0) {
                classes <- c('kingdom','superclass','class','subclass') %>%
                  {.[. %in% names(classi)]}
                
                classi <- classi %>%
                  select(INCHIKEY,{{classes}},contains('level')) %>%
                  left_join(x %>%
                              hits() %>%
                              entries() %>%
                              select(ID,INCHIKEY), by = "INCHIKEY") %>%
                  select(ID,INCHIKEY,everything())
                
                message(str_c(length(unique(classi$INCHIKEY)),' classifications returned'))  
              } else 
                message('0 classifications returned')
              
              x@classifications <- classi  
            }
            
            return(x)
          }
)
