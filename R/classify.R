
#' @importFrom magrittr set_names
#' @importFrom dplyr distinct select
#' @importFrom classyfireR get_classification
#' @importFrom purrr map_lgl
#' @importFrom tidyr spread

classify <- function(inchikey){
  
  suppressMessages(
    classi <- inchikey %>%
      map(get_classification) %>%
      set_names(inchikey)
  )
  
  classi %>%
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
}