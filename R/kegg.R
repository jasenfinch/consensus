#' keggConsensus
#' @rdname keggConsensus
#' @description Collate consensus classifications for molecular formula assignments using KEGG.
#' @param x S4 object of class Assignment
#' @param organism organism kegg ID.
#' @param threshold majority assignment threshold for consensus classifications
#' @importFrom FELLA buildGraphFromKEGGREST
#' @importFrom tidygraph as_tbl_graph
#' @importFrom MFassign nodes
#' @importFrom dplyr anti_join full_join rowwise
#' @importFrom utils capture.output
#' @importFrom mzAnnotation descriptors
#' @export

setMethod('keggConsensus',signature = 'Assignment',
          function(x,organism = 'hsa', threshold = 0.5){
            
            adductRules <- x@parameters@adductRules
            
            capture.output({
              suppressMessages({
                g <- buildGraphFromKEGGREST(organism = organism, filter.path = NULL) %>%
                  as_tbl_graph() %>%
                  nodes() %>%
                  filter(com == 5)
              })
            })
            
            suppressMessages(met <- metabolites %>%
              filter(ACCESSION_ID %in% g$name) %>%
              {metaboliteDB(.,descriptors = descriptors(.))})
            
            mfs <- x %>%
              assignments() %>%
              select(Name,MF,Adduct) %>%
              distinct()
            
            hits <- mfs %>%
              split(1:nrow(.)) %>%
              map(~{
                m <- .
                met %>%
                  filterMF(m$MF) %>%
                  getAccessions() %>%
                  mutate(Name = m$Name,MF = m$MF,Adduct = m$Adduct)
              }) %>%
              bind_rows() %>%
              rowwise() %>%
              mutate(inchikey = mzAnnotation::convert(INCHI,'inchi','inchikey'))
            
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
              mutate(inchikey = mzAnnotation::convert(INCHI,'inchi','inchikey'))
            
            noPIPs <- mfs %>%
              select(-Name) %>%
              anti_join(pips, by = c("MF", "Adduct")) %>%
              mutate(kingdom = 'No hits')
            
            classifications <- pips$inchikey %>%
              unique() %>%
              classify() %>%
              right_join(pips,by = c('InChIKey' = 'inchikey')) %>%
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