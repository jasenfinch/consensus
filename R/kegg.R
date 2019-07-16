#' @importFrom FELLA buildGraphFromKEGGREST
#' @importFrom tidygraph as_tbl_graph
#' @importFrom MFassign nodes
#' @importFrom dplyr anti_join full_join rowwise
#' @importFrom utils capture.output
#' @importFrom mzAnnotation descriptors
#' @export

setMethod('keggConsensus',signature = 'Assignment',
          function(x,organism = 'hsa', threshold = 0.5){
            
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
            
            MFhits <- mfs %>%
              split(1:nrow(.)) %>%
              map(~{
                m <- .
                met %>%
                  filterMF(m$MF) %>%
                  filterIP(Adducts$Rule[Adducts$Name == m$Adduct]) %>%
                  getAccessions() %>%
                  mutate(Name = m$Name,MF = m$MF,Adduct = m$Adduct)
              }) %>%
              bind_rows() %>%
              rowwise() %>%
              mutate(inchikey = mzAnnotation::convert(INCHI,'inchi','inchikey'))
            
            noHits <- mfs %>%
              select(-Name) %>%
              anti_join(MFhits, by = c("MF", "Adduct")) %>%
              mutate(kingdom = 'No hits')
            
            classifications <- MFhits$inchikey %>%
              unique() %>%
              classify() %>%
              right_join(MFhits,by = c('InChIKey' = 'inchikey')) %>%
              select(ACCESSION_ID:last_col(),everything())
            
            noClassifications <- MFhits %>%
              anti_join(classifications, by = c("ACCESSION_ID", "MF", "INCHI", "SMILE", "Name", "Adduct")) %>%
              select(MF,Adduct) %>%
              mutate(kingdom = 'unclassified')
            
            consensi <- classifications %>%
              consensusCls(threshold = threshold) %>%
              full_join(noHits, by = c("MF", "Adduct", "kingdom")) %>%
              full_join(noClassifications, by = c("MF", "Adduct", "kingdom"))
            return(consensi)
          }
)