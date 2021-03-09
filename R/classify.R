
#' @importFrom magrittr set_names %>%
#' @importFrom dplyr distinct select bind_rows filter left_join
#' @importFrom classyfireR get_classification classification
#' @importFrom purrr map_dbl map
#' @importFrom tidyr spread
#' @importFrom tidyselect last_col
#' @importFrom progress progress_bar
#' @importFrom stringr str_c
#' @importFrom tibble is_tibble

setMethod('classify',signature = 'Consensus',
          function(x){
            inchikey <- x %>%
              hits() %>%
              getAccessions() %>%
              .$INCHIKEY
            
            if (length(inchikey) > 0) {
              message(str_c('Retrieving classifications for ',length(inchikey),' InChIKeys...'))
              
              pb <- progress_bar$new(
                format = "[:bar] :percent eta: :eta",
                total = length(inchikey), clear = FALSE)
              pb$tick(0)
              
              classi <- inchikey %>%
                map(~{
                  cl <- suppressMessages(get_classification(.)) 
                  if (is.null(cl)) {
                    cl <- tibble(Level = 'kingdom','Classification' = 'Unclassified',CHEMONT = NA) 
                  } else {
                    cl <- cl %>%
                      classification()
                  }
                  Sys.sleep(5)
                  pb$tick()
                  return(cl)
                }) %>%
                set_names(inchikey)
              
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
                              getAccessions() %>%
                              select(ID,INCHIKEY), by = "INCHIKEY") %>%
                  select(ID,INCHIKEY,everything())
                
                message(str_c(length(unique(classi$INCHIKEY)),' classifications returned'))  
              } else {
                message('0 classifications returned')
              }
              
              x@classifications <- classi  
            }
            
            return(x)
          }
)