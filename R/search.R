#' @importFrom mzAnnotation filterMF metaboliteDB filterACCESSIONS descriptors convert
#' @importFrom dplyr rowwise ungroup

setMethod('mfHits',signature = 'Consensus',
          function(x){
            
            db <- database(x)
            
            if (db == 'kegg') {
              message('Searching KEGG...')
              compounds <- keggCompounds(organism(x))
              
              met <- metabolites %>%
                filterACCESSIONS(compounds)
              
              hits <- met %>%
                filterMF(mf(x)) %>%
                {
                  if (nrow(getAccessions(.)) > 0) {
                    .@accessions[[1]] <- .@accessions[[1]] %>%
                      rowwise() %>%
                      mutate(INCHIKEY = convert(INCHI,'inchi','inchikey')) %>%
                      ungroup() 
                  } else {
                    .@accessions[[1]] <- .@accessions[[1]] %>%
                      mutate(INCHIKEY = character())
                  }
                  
                  return(.)
                }
              
              message(str_c(hits %>% 
                              getAccessions() %>% 
                              nrow()), ' hits returned')
            }
            
            if (db == 'pubchem') {
              hits <- pubchemMatch(mf(x)) %>%
                {metaboliteDB(.,descriptors(.))}
            }
            
            x@hits <- hits
            return(x)
          })

#' @importFrom mzAnnotation filterIP

setMethod('pips',signature = 'Consensus',
          function(x){
            
            a <- adductRules(x)
            h <- hits(x)
            
            p <- a$Rule %>%
              map(filterIP,db = h) %>%
              set_names(a$Name) %>%
              map(getAccessions) %>%
              bind_rows(.id = 'Adduct') %>%
              select(Adduct,ACCESSION_ID)
            
            
            x@PIPs <- p
            
            return(x)
          })