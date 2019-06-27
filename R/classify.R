
#' @importFrom magrittr set_names
#' @importFrom dplyr distinct select
#' @importFrom classyfireR get_classification
#' @importFrom purrr map_chr
#' @importFrom tidyr spread
#' @importFrom parallel detectCores makeCluster stopCluster parLapply

classify <- function(inchikey,nCores = availableCores() * 0.75){
  
  # plan(multiprocess,workers = nCores)
  # 
  # suppressMessages(classi <- inchikey %>%
  #   future_map(get_classification) %>%
  #   set_names(inchikey))
  
  # clus <- makeCluster(nCores,type = 'FORK')
  # suppressMessages(classi <- inchikey %>%
  #                    parLapply(cl = clus,fun = get_classification) %>%
  #                    set_names(inchikey)
  # )
  # stopCluster(clus)
  
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