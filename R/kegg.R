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
            select(ENTRY,NAME,FORMULA,EXACT_MASS,PATHWAY,REACTION,ENZYME) %>%
            rename(SYNONYM = NAME) %>%
            mutate(NAME = str_split_fixed(SYNONYM,';',2)[,1],
                   SYNONYM = str_split_fixed(SYNONYM,';',2)[,2]) %>%
              select(ENTRY,NAME,SYNONYM:ENZYME)
          }) %>%
        bind_rows()
    }) %>%
    bind_rows()
}

#' keggCompounds
#' @description Return KEGG compound accession IDs for a given organism
#' @param organism KEGG organism ID. If NULL use all compounds.
#' @examples 
#' \dontrun{
#' keggCompounds('hsa')
#' }
#' @importFrom KEGGREST keggLink keggGet
#' @importFrom stringr str_remove_all coll str_split_fixed str_split
#' @importFrom tibble deframe
#' @export

keggCompounds <- function(organism = NULL){
  
  if (is.null(organism)) {
    compounds <- metabolites$ACCESSION_ID
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

#' keggPIPs
#' @description Return KEGG putative ionisation products (PIPs) for given molecular formulas and their adducts.
#' @param MFs tibble containing 2 columns. The first named MF containing molecular formulas. The second named Adduct containing adducts. See example.
#' @param organism KEGG organism ID. If NULL use all compounds.
#' @examples 
#' \dontrun{
#' d <- tibble(MF = c('C12H22O11','C4H6O5'),Adduct = rep('[M-H]1-',2))
#' keggPIPs(d,'hsa')
#' }
#' @importFrom dplyr tbl_df
#' @importFrom mzAnnotation filterACCESSIONS adducts convert
#' @export

keggPIPs <- function(MFs,organism = NULL, adductRules = adducts()){
  
  compounds <- keggCompounds(organism)
  
  met <- metabolites %>%
    filterACCESSIONS(compounds)
  
  hits <- MFs %>%
    split(1:nrow(.)) %>%
    map(~{
      m <- .
      met %>%
        filterMF(m$MF) %>%
        getAccessions() %>%
        mutate(MF = m$MF,Adduct = m$Adduct)
    }) %>%
    bind_rows() %>%
    mutate(Name = str_c(MF,' ',Adduct)) %>%
    rowwise() %>%
    select(Name,MF,Adduct,everything())
  
  pips <- hits %>%
    split(.$Name) %>%
    map(~{
      m <- .
      met %>%
        filterMF(m$MF[1]) %>%
        filterIP(adductRules$Rule[adductRules$Name == m$Adduct[1]]) %>%
        getAccessions() %>%
        mutate(Name = m$Name[1],MF = m$MF[1],Adduct = m$Adduct[1])
    }) %>%
    bind_rows() %>%
    rowwise() %>%
    mutate(INCHIKEY = convert(INCHI,'inchi','inchikey')) %>%
    tbl_df() %>%
    select(MF,Adduct,everything(),-Name)
  
  return(pips)
}

#' keggConsensus
#' @rdname keggConsensus
#' @description Collate consensus classifications for molecular formula assignments using KEGG.
#' @param x S4 object of class Assignment
#' @param organism organism kegg ID. If NULL use all compounds.
#' @param threshold majority assignment threshold for consensus classifications
#' @importFrom dplyr anti_join full_join rowwise
#' @importFrom mzAnnotation descriptors
#' @export

setMethod('keggConsensus',signature = 'Assignment',
          function(x,organism = NULL, threshold = 0.5){
            
            adductRules <- x@parameters@adductRules
            
            mfs <- x %>%
              assignments() %>%
              select(Name,MF,Adduct) %>%
              distinct()
            
            pips <- mfs %>%
              select(MF,Adduct) %>%
              distinct() %>%
              keggPIPs()
            
            noPIPs <- mfs %>%
              select(-Name) %>%
              anti_join(pips, by = c("MF", "Adduct")) %>%
              mutate(kingdom = 'No hits')
            
            classifications <- pips$INCHIKEY %>%
              unique() %>%
              classify() %>%
              right_join(pips,by = c('InChIKey' = 'INCHIKEY')) %>%
              select(ACCESSION_ID:last_col(),everything())
            
            noClassifications <- pips %>%
              anti_join(classifications, by = c("ACCESSION_ID", "MF", "INCHI", "SMILE", "Name", "Adduct")) %>%
              select(MF,Adduct) %>%
              mutate(kingdom = 'Unclassified')
            
            consensi <- new('Consensus')
            consensi@hits <- hits
            consensi@PIPs <- pips
            consensi@classifications <- classifications
            consensi@consensus <- classifications %>%
              consensusCls(threshold = threshold) %>%
              full_join(noPIPs, by = c("MF", "Adduct", "kingdom")) %>%
              full_join(noClassifications, by = c("MF", "Adduct", "kingdom"))
            
            return(consensi)
          }
)