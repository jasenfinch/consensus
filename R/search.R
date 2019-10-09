#' @importFrom mzAnnotation filterMF

setMethod('mfHits',signature = 'Consensus',
          function(x){
            db <- database(x)
            
            if (db == 'kegg') {
              message('Searching KEGG...')
              compounds <- keggCompounds(organism(x))
              
              met <- metabolites %>%
                filterACCESSIONS(compounds)
              
              hits <- met %>%
                filterMF(mf(x))
              message(str_c(hits %>% 
                        getAccessions() %>% 
                        nrow()), 'hits returned')
            }
            
            if (db == 'pubchem') {
              hits <- pubchemMatch(MF) %>%
                {metaboliteDB(.,descriptors(.))}
            }
            
            x@hits <- hits
            return(x)
          })