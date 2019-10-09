#' pubchemMatch
#' @description Query the pubchem database for molecular formula matches.
#' @param MF a molecular formula to query
#' @importFrom dplyr rename mutate
#' @importFrom httr GET content config
#' @importFrom tibble as_tibble
#' @examples
#' \dontrun{ 
#' pubchemMatch('C12H22O11')
#' }
#' @export

pubchemMatch <- function(MF){
  message(str_c('\n',MF))
  message('Retreiving CIDs...')
  cmd <- str_c('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/formula/',
               MF,
               '/JSON')
  
  key <- cmd %>%
    GET(config(http_version = 0)) %>%
    content() %>%
    {.$Waiting$ListKey}
  
  while (T) {
    Sys.sleep(0.3)
    cid_cmd <- str_c('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/',
                     key,
                     '/cids/JSON')
    
    cids <- cid_cmd %>%
      GET(config(http_version = 0)) %>%
      content()
    
    if (names(cids) == 'Waiting') {
      key <- cids$Waiting$ListKey
    } else {
      break()
    }
  }
  
  if (names(cids) == 'Fault') {
    message('0 CIDs returned')
    return(NULL)
  }
  
  if (names(cids) == "IdentifierList") {
    cids <- cids %>%
      {.$IdentifierList$CID} %>%
      unlist()
    
    if (length(cids) > 200) {
      cids <- cids %>%
        split(., ceiling(seq_along(.)/200))
    } else {
      cids <- list(cids)
    }
    
    descs <- c('IUPACName','MolecularFormula','Charge','InChI','InChIKey','CanonicalSMILES','CovalentUnitCount')
    chem_info_cmd <- cids %>%
      map(~{
        str_c('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',
              str_c(.,collapse = ','),
              '/property/',
              str_c(descs,collapse = ','),
              '/JSON')
      })
    
    chem_info <- chem_info_cmd %>%
      map(~{
        Sys.sleep(0.3)
        GET(.,config(http_version = 0)) %>%
          content() %>%
          {.$PropertyTable$Properties} %>%
          map(as_tibble) %>%
          bind_rows()
      })  %>%
      bind_rows() %>%
      filter(Charge == 0)
    
    if (nrow(chem_info) > 0) {
     chem_info <- chem_info  %>%
        mutate(ACCESSION_ID = 1:nrow(.)) %>%
        rename(SMILE = CanonicalSMILES) %>%
        filter(CovalentUnitCount == 1) %>%
       rename(MF = MolecularFormula)
    } else {
      message('0 CIDs returned')
      return(NULL)
    }
    
   
  }
  message(str_c(nrow(chem_info),' CIDs returned'))
  return(chem_info)
}

#' pubchemConsensus
#' @description Calculate a consensus classification for a given molecular formula and adducts based on the pubchem database.
#' @param MF molecular formula
#' @param adducts character vector of adducts
#' @param threshold consensus threshold
#' @param adductRules Adduct formation rules to use for putative ionisation products. Defaults to \code{mzAnnotation::adducts()}
#' @examples
#' \dontrun{
#' pubchemConsensus('C10H10O7')
#' }
#' @importFrom stringr str_detect
#' @importFrom tibble tibble
#' @export

pubchemConsensus <- function(MF, adducts = c('[M-H]1-'), threshold = 0.5, adductRules = adducts()){
  hits <- pubchemMatch(MF)
  
  if (is.null(hits)) {
    hits <- tibble()
    PIPs <- tibble()
    classifications <- tibble()
    con <-  tibble(MF = MF,Adduct = adducts,kingdom = 'No hits')
  } else {
    PIPs <- pips(hits,adducts)
    
    if (nrow(PIPs) == 0) {
      PIPs <- tibble()
      classifications <- tibble()
      con <- tibble(MF = MF,Adduct = adducts,kingdom = 'No hits')
    } else {
      classifications <- pipClassifications(PIPs)
      
      if (nrow(classifications) == 0) {
        classifications <- tibble()
        con <- tibble(MF = MF,Adduct = adducts,kingdom = 'Unclassified')
      } else {
        con <- classifications %>%
          consensus(threshold = threshold)   
      }
      
    }
  }
  
  consensus <- new('Consensus')
  consensus@hits <- hits
  consensus@PIPs <- PIPs
  consensus@classifications <- classifications
  consensus@consensus <- con
  
  return(consensus)
}