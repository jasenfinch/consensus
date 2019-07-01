
#' @importFrom magrittr set_names
#' @importFrom dplyr distinct select
#' @importFrom classyfireR get_classification
#' @importFrom purrr map_lgl
#' @importFrom tidyr spread
#' @importFrom tidyselect last_col
#' @import progress

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
    filter(!is.na(kingdom)) %>%
    select(InChIKey,kingdom,superclass,class,subclass,`level 5`:last_col())
  
  message(str_c(length(unique(classi$InChIKey)),' classifications returned'))
  
  return(classi)
}