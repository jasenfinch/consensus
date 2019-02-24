
#' @importFrom magrittr set_names
#' @importFrom dplyr distinct select
#' @importFrom classyfireR entity_classification
#' @importFrom purrr map_chr
#' @importFrom tidyr spread


classify <- function(inchis){
  classi <- inchis %>%
    map_chr(convert,inputType = 'inchi',outputType = 'inchikey') %>%
    map(entity_classification) %>%
    set_names(inchis)
  
  # missing <- names(classifications)[sapply(classifications,is.null)]
  
  classi %>%
    .[!sapply(.,is.null)] %>%
    map(~{
      d <- .
      d %>%
        select(-CHEMONT) %>%
        spread(.,Level,Classification)
    }) %>%
    bind_rows(.id = 'InChI') %>%
    distinct() %>%
    filter(!is.na(kingdom)) %>%
    select(InChI,kingdom,superclass,class,subclass)
  
}