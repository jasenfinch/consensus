
#' @importFrom magrittr set_names %>%
#' @importFrom dplyr distinct select bind_rows filter left_join
#' @importFrom classyfireR get_classification
#' @importFrom purrr map_lgl map
#' @importFrom tidyr spread
#' @importFrom tidyselect last_col
#' @importFrom progress progress_bar
#' @importFrom stringr str_c

setMethod('classify',signature = 'Consensus',
          function(x){
            inchikey <- x %>%
              hits() %>%
              getAccessions() %>%
              .$INCHIKEY
            
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
                }
                Sys.sleep(5)
                pb$tick()
                return(cl)
              }) %>%
              set_names(inchikey)
            
            classes <- c('kingdom','superclass','class','subclass')
            
            classi <- classi %>%
              .[!map_lgl(.,is.null)] %>%
              map(~{
                d <- .
                d %>%
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
                            select(ACCESSION_ID,INCHIKEY), by = "INCHIKEY") %>%
                select(ACCESSION_ID,INCHIKEY,everything())
              
              message(str_c(length(unique(classi$INCHIKEY)),' classifications returned'))  
            } else {
              message('0 classifications returned')
            }
            
            x@classifications <- classi
            
            return(x)
          }
)