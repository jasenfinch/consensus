#' @importFrom FELLA buildGraphFromKEGGREST
#' @importFrom tidygraph as_tbl_graph
#' @importFrom MFassign nodes

setMethod('keggConsensus',signature = 'Assignment',
          function(x,organism = 'hsa'){
            g <- suppressMessages(buildGraphFromKEGGREST(organism = organism, filter.path = NULL) %>%
                as_tbl_graph() %>%
                nodes() %>%
                filter(com == 5))
            
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
            
            classifications <- MFhits$inchikey %>%
              unique() %>%
              classify() %>%
              left_join(MFhits,by = c('InChIKey' = 'inchikey')) %>%
              rename(MolecularFormula = MF)
            
            classifications %>%
              consensusCls(threshold = threshold) %>%
              select('MolecularFormula','Adduct','Score',everything())
          }
)