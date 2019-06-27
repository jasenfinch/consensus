

pubchemMatch <- function(MF){
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
        filter(CovalentUnitCount == 1)  
    } else {
      message('0 CIDs returned')
      return(NULL)
    }
    
   
  }
  message(str_c(nrow(chem_info),' CIDs returned'))
  return(chem_info)
}

pips <- function(hits,adducts){
  descriptorTable <- descriptors(hits)
  db <- metaboliteDB(hits,descriptorTable)
  
  accessions <- adducts %>%
    map(~{
      rule <- Adducts$Rule[Adducts$Name == .]
      
      ips <- filterIP(db,rule)
      
      ips <- ips %>%
        getAccessions()
      message(nrow(ips),' of ',nrow(hits),' can ionize as ',.)
      return(ips)
    }) %>%
    set_names(adducts) %>%
    bind_rows(.id = 'Adduct')
  
  return(accessions)
}


#' pubchemPIPs
#' @examples 
#' pubchemPIPs('C4H6O5', '[M-H]1-')
#' @importFrom stringr str_c
#' @importFrom httr GET content config
#' @importFrom purrr map
#' @importFrom dplyr bind_rows filter mutate rowwise
#' @importFrom mzAnnotation descriptors metaboliteDB getAccessions 
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @importFrom progress progress_bar

setMethod('pubchemPIPs',signature = 'Consensus',
          function(ips){
            iprods <- ips@IPs
            
            pb <- progress_bar$new(
              format = "  retrieving [:bar] :percent eta: :eta",
              total = length(unique(iprods$MF)), clear = FALSE)
            pb$tick(0)
            
            PIPs <- iprods %>%
              split(.$MF) %>%
              map(~{
                d <- .
                MF <- d$MF[1]
                IPs <- d$Adduct
                
                message(str_c('\n',MF))
                
                hits <- pubchemMatch(MF)
                
                if (!is.null(hits)) {
                  accessions <- pips(hits,IPs)
                } else {
                  return(NULL)
                }
                
                if (!is.null(accessions)) {
                  accessions <- accessions %>%
                    select(CID,MolecularFormula,Adduct,SMILE:Charge)
                } else {
                  return(NULL)
                }
                pb$tick()
                message('\n')
                return(accessions)
              }) %>%
              bind_rows()
            ips@PIPs <- PIPs
            return(ips)
          }
)