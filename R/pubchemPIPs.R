#' pubchemPIPs
#' @example 
#' pubchemPIPs('C4H6O5', '[M-H]1-')
#' @importFrom stringr str_c
#' @importFrom httr GET content config
#' @importFrom purrr map
#' @importFrom dplyr bind_rows filter mutate rowwise
#' @importFrom mzAnnotation descriptors metaboliteDB getAccessions 
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble

pubchemPIPs <- function(MF,IP){
  message(str_c('\n',MF,' ',IP,'\n'))
  
  cmd <- str_c('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/formula/',
               MF,
               '/JSON')
  
  key <- cmd %>%
    GET(config(http_version = 0)) %>%
    content() %>%
    {.$Waiting$ListKey}
  
  while (T) {
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

    message(str_c(length(cids),' CIDs returned...'))
    
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
      filter(Charge == 0) %>%
      mutate(ACCESSION_ID = 1:nrow(.)) %>%
      rename(SMILE = CanonicalSMILES) %>%
      filter(CovalentUnitCount == 1)
    
    descriptorTable <- mzAnnotation::descriptors(chem_info)
    
    db <- metaboliteDB(chem_info,descriptorTable)
    
    rule <- Adducts$Rule[Adducts$Name == IP]
    
    ips <- mzAnnotation:::filterIP(db,rule)
    
    accessions <- ips %>%
      mzAnnotation::getAccessions()
    
    message('...of which ',nrow(accessions),' can ionize')
    
    return(accessions) 
  }
}