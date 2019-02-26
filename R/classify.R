
#' @importFrom magrittr set_names
#' @importFrom dplyr distinct select
#' @importFrom classyfireR get_classification
#' @importFrom purrr map_chr
#' @importFrom tidyr spread


classify <- function(inchikey){
  classi <- inchikey %>%
    map(get_classification) %>%
    set_names(inchikey)
  
  # missing <- names(classifications)[sapply(classifications,is.null)]
  
  classi %>%
    .[!sapply(.,is.null)] %>%
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