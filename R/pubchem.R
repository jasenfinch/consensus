#' @importFrom httr GET config content
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate rename

pubchemMatch <- function(MF){
  message('Retrieving CIDs from PubChem...')
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
    return(tibble(ID = character(),
                  NAME = character(),
                  MF = character(),
                  INCHI = character(),
                  SMILES = character(),
                  INCHIKEY = character()))
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
        mutate(ID = 1:nrow(.)) %>%
        rename(NAME = IUPACName,
               INCHI = InChI,
               SMILES = CanonicalSMILES, 
               INCHIKEY = InChIKey) %>%
        filter(CovalentUnitCount == 1) %>%
       rename(MF = MolecularFormula)
    } else {
      message('0 CIDs returned')
      return(tibble(CID = integer(),
             MF = character(),
             SMILES = character(),
             INCHI = character(),
             INCHIKEY = character(),
             NAME = character(),
             Charge = integer(),
             CovalentUnitCount = integer(),
             ID = integer()))
    }
    
   
  }
  message(str_c(nrow(chem_info),' CIDs returned'))
  return(chem_info)
}
