#' @importFrom mzAnnotation filterMF metaboliteDB filterACCESSIONS descriptors convert
#' @importFrom dplyr rowwise tbl_df

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
                              nrow()), ' hits returned')
            }
            
            if (db == 'pubchem') {
              hits <- pubchemMatch(MF) %>%
                {metaboliteDB(.,descriptors(.))}
            }
            
            if (nrow(hits@accessions[[1]]) > 0) {
              hits@accessions[[1]] <- hits@accessions[[1]] %>%
                rowwise() %>%
                mutate(INCHIKEY = convert(INCHI,'inchi','inchikey')) %>%
                tbl_df() 
            } else {
              hits@accessions[[1]] <-  hits@accessions[[1]] %>%
                mutate(INCHIKEY = character())
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