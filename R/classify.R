
#' @importFrom magrittr set_names %>%
#' @importFrom dplyr distinct select bind_rows filter
#' @importFrom classyfireR get_classification
#' @importFrom purrr map_lgl map
#' @importFrom tidyr spread
#' @importFrom tidyselect last_col
#' @importFrom progress progress_bar
#' @importFrom stringr str_c

classify <- function(inchikey){
  
  message(str_c('Retreiving classifications for ',length(inchikey),' InChIKeys...'))
  
  pb <- progress_bar$new(
    format = "[:bar] :percent eta: :eta",
    total = length(inchikey), clear = FALSE)
  pb$tick(0)
  
  classi <- inchikey %>%
    map(~{
      cl <- suppressMessages(get_classification(.))
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
    bind_rows(.id = 'InChIKey') %>%
    distinct() %>%
    filter(!is.na(kingdom))
  
  if (nrow(classi) > 0) {
    classes <- c('kingdom','superclass','class','subclass') %>%
      {.[. %in% names(classi)]}
    
    classi <- classi %>%
      select(InChIKey,{{classes}},contains('level'))
    
    message(str_c(length(unique(classi$InChIKey)),' classifications returned'))  
  } else {
    message('0 classifications returned')
  }
  
  return(classi)
}