
#' @importFrom magrittr set_names
#' @importFrom dplyr distinct select
#' @importFrom classyfireR get_classification
#' @importFrom purrr map_chr
#' @importFrom tidyr spread
#' @importFrom parallel detectCores makeCluster stopCluster parLapply
#' @export

classify <- function(inchikey,nCores = detectCores(), clusterType = 'PSOCK'){
  clus <- makeCluster(nCores,type = clusterType)
  classi <- inchikey %>%
    parLapply(cl = clus,fun = get_classification) %>%
    set_names(inchikey)
  stopCluster(clus)
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