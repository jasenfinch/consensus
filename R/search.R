#' @importFrom mzAnnotation filterMF metaboliteDB filterEntries descriptors convert
#' @importFrom dplyr rowwise ungroup

setMethod('mfHits',signature = 'Consensus',
          function(x){
            
            db <- database(x)
            
            if (db == 'kegg') {
              message('Searching KEGG...')
              compounds <- keggCompounds(organism(x))
              
              met <- metabolites %>%
                filterEntries(compounds)
              
              hits <- met %>%
                filterMF(mf(x)) %>%
                {
                  if (nrow(entries(.)) > 0) {
                    entries(.) <- entries(.) %>%
                      rowwise() %>%
                      mutate(INCHIKEY = convert(INCHI,'inchi','inchikey')) %>%
                      ungroup() 
                  } else {
                    entries(.) <- entries(.) %>%
                      mutate(INCHIKEY = character())
                  }
                  
                  .
                }
              
              message(str_c(hits %>% 
                              entries() %>% 
                              nrow()), ' hits returned')
            }
            
            if (db == 'pubchem') {
              hits <- pubchemMatch(mf(x)) %>%
                {metaboliteDB(.,descriptors(.$SMILES))}
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
              map(entries) %>%
              bind_rows(.id = 'Adduct') %>%
              select(Adduct,ID)
            
            
            x@PIPs <- p
            
            return(x)
          })