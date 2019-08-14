globalVariables(c('.','kingdom','CID','MF','Adduct','InChIKey','superclass','subclass','level 5','MF',
                  'Charge','CanonicalSMILES','CovalentUnitCount','SMILE','ACCESSION_ID','MolecularFormula','com','Name','INCHI','Score','Feature','Intensity'))

#' construct
#' @rdname construct
#' @description Consensus structural classifications for molecular formula assignments.
#' @param x S4 object of class Workflow
#' @param organism organism kegg ID.
#' @param threshold majority assignment threshold for consensus classifications
#' @param databases databases to use for PIP searches. Either kegg or pubchem.
#' @importClassesFrom MFassign Assignment
#' @importFrom MFassign assignments
#' @importFrom methods new
#' @importFrom lubridate seconds_to_period
#' @examples
#' \dontrun{
#' library(MFassign) 
#' p <- assignmentParameters('FIE')
#' p@nCores <- 2
#' assignment <- assignMFs(peakData,p)
#' 
#' consensusCl <- construct(assignment)
#' }
#' @export

setMethod('construct',signature = 'Assignment',
          function(x, organism = 'hsa', threshold = 0.5, databases = c('kegg','pubchem')){
            
            consense <- new('Consensuses')
            
            adductRules <- x@parameters@adductRules
            
            if ('kegg' %in% databases) {
              z <- keggConsensus(x,organism = organism)
              
              n <- z %>%
                .@consensus %>%
                filter(kingdom == 'No hits' | kingdom == 'Unclassified') %>%
                distinct() 
            } else {
              
            }
            
            if ('pubchem' %in% databases | nrow(n) > 0) {
              startTime <- proc.time()
              pc <- n %>%
                select(MF,Adduct) %>%
                split(.$MF) %>%
                map(~{
                  consensusClassification(.$MF[1],.$Adduct,adductRules = adductRules)
                }) 
              endTime <- proc.time()
              
              elapsed <- {endTime - startTime} %>%
                .[3] %>%
                round(1) %>%
                seconds_to_period() %>%
                str_c('[',.,']')
              
              cat('\n',elapsed)
              con <- pc %>%
                map(~{
                  .@consensus
                }) %>%
                bind_rows() %>%
                select(MF,Adduct,Score,everything()) %>%
                bind_rows(z@consensus %>%
                            filter(kingdom != 'No hits' & kingdom != 'Unclassified'))
              
              dat <- x %>%
                .@data %>%
                gather('Feature','Intensity') %>%
                group_by(Feature) %>%
                summarise(Intensity = mean(Intensity))
              dat[dat == ''] <- NA
              dat <- dat %>%
                left_join(con %>%
                            left_join(x %>% 
                                        assignments() %>% 
                                        select(Name,Feature,MF,Adduct), 
                                      by = c("MF", "Adduct")), 
                          by = c('Feature')) %>%
                select(Name,everything())
              dat$Name[is.na(dat$Name)] <- dat$Feature[is.na(dat$Name)]
              dat <- dat %>%
                select(-Feature)
              dat$kingdom[is.na(dat$kingdom)] <- 'Unknown'
              
              consense@consensuses <- c(list(KEGG = z),pc)
              consense@results <- dat 
            } else {
              consense@consensuses <- list(KEGG = z)
              consense@results <- z@consensus
            }
            
            return(consense)
          }
)