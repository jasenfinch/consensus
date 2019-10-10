#' keggCompoundInfo
#' @description Return KEGG compound accession information.
#' @param IDs KEGG compound accession IDs
#' @examples 
#' keggCompoundInfo(c('C00089','C00149'))
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
            as_tibble() %>%
            select(ENTRY,NAME,FORMULA,EXACT_MASS) %>%
            rename(SYNONYM = NAME) %>%
            mutate(NAME = str_split_fixed(SYNONYM,';',2)[,1],
                   SYNONYM = str_split_fixed(SYNONYM,';;',2)[,2]) %>%
              select(ENTRY,NAME,SYNONYM:EXACT_MASS)
          }) %>%
        bind_rows()
    }) %>%
    bind_rows()
}

#' @importFrom KEGGREST keggLink
#' @importFrom tibble deframe tibble
#' @importFrom stringr str_remove_all

keggCompounds <- function(organism = character()){
  
  if (length(organism) == 0) {
    compounds <- metabolites %>%
      getAccessions() %>%
      .$ACCESSION_ID
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
