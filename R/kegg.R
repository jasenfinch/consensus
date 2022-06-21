#' keggCompoundInfo
#' @description Return KEGG compound accession information.
#' @param IDs KEGG compound accession IDs
#' @examples 
#' keggCompoundInfo(c('C00089','C00149'))
#' @importFrom KEGGREST keggGet
#' @importFrom stringr str_split_fixed
#' @export

keggCompoundInfo <- function(IDs){
  IDs %>%
    split(ceiling(seq_along(.)/10)) %>%
    map(~{
      keggGet(.) %>%
        map(~{
          d <- .
          d %>%
            map(~{str_c(.,collapse = '; ')}) %>%
            as_tibble()
          }) %>%
        bind_rows()
    }) %>%
    bind_rows()
}

#' @importFrom KEGGREST keggLink
#' @importFrom tibble deframe tibble
#' @importFrom stringr str_remove_all
#' @importFrom mzAnnotation entries

keggCompounds <- function(organism = character()){
  
  if (length(organism) == 0) {
    compounds <- metabolites %>%
      entries() %>%
      .$ID
  } else {
    enzymes <- keggLink(organism,'enzyme') %>%
      names()
    compounds <- keggLink('compound','enzyme') %>%
      {tibble(Enzyme = names(.),Compound = .)} %>%
      filter(Enzyme %in% enzymes) %>%
      select(Compound) %>%
      distinct() %>%
      deframe() %>%
      str_remove_all('cpd:')
  }
  return(compounds)
}
