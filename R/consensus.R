
globalVariables(c('.','kingdom','CID','MF','Adduct','InChIKey','superclass','subclass','level 5','MF',
                  'Charge','CanonicalSMILES','CovalentUnitCount','SMILE','ACCESSION_ID','MolecularFormula','com','Name','INCHI'))

#' @importFrom tibble rowid_to_column
#' @importFrom dplyr everything group_by summarise right_join
#' @importFrom tidyr gather
#' @importFrom tidyselect contains

consensusCls <- function(classifications,threshold = 0.5){
  
  if (nrow(classifications) > 1) {
    consensusClasses <- classifications %>%
      split(str_c(.$MF,.$Adduct,sep = ' ')) %>%
      map(~{
        cl <- .
        levels <- names(cl)[which(names(cl) == 'kingdom'):length(names(cl))]
        
        classes <- cl %>%
          select(MF,Adduct,everything()) %>%
          distinct() %>%
          rowid_to_column(var = 'ID')
        
        suppressMessages(freq <- cl %>%
                           left_join(classes) %>%
                           group_by(ID) %>%
                           summarise(N = n()) %>%
                           right_join(classes, by = "ID"))
        
        votes <- levels %>%
          map(~{
            lev <- .
            freq %>%
              rename('Class' = !!lev) %>%
              select(ID:Adduct,Class) %>%
              split(.$Class) %>%
              map(~{
                d <- .
                d %>%
                  group_by(MF,Adduct,Class) %>%
                  summarise(N = sum(N))
              }) %>%
              bind_rows()
          }) %>%
          set_names(levels) %>%
          bind_rows(.id = 'Level')
        
        clLevels <- c('kingdom','superclass','class','subclass')
        
        votesTable <- freq %>%
          select(-N) %>%
          split(1:nrow(.)) %>%
          map(~{
            d <- .
            d %>%
              gather('Level','Class',-(ID:Adduct)) %>%
              left_join(votes, by = c("MF", "Adduct", "Level", "Class")) %>%
              select(-Class) %>%
              spread(Level,N)
          }) %>%
          bind_rows() 
        
        clLevels <- clLevels[clLevels %in% names(votesTable)]
        
        votesTable <- votesTable %>%
          select(ID:Adduct,clLevels,contains('level'))
        
        N <- nrow(cl)
        
        p <- votesTable %>%
          mutate(N = N) %>%
          select(N,kingdom:names(.)[length(names(.))])
        proportions <- p
        for (i in 2:ncol(proportions)) {
          proportions[,i] <- p[,i] / p[,i - 1]
        }
        proportions <- proportions %>%
          select(-N) %>%
          mutate(ID = 1:nrow(.))
        
        consensus <- proportions %>%
          select(-ID) %>% 
          split(1:nrow(.)) %>%
          map(~{
            mutate(.,Score = prod(.,na.rm = T))  
          }) %>%
          bind_rows() %>%
          mutate(ID = 1:nrow(.))
        
        maxScore <- max(consensus$Score)
        
        cons <- consensus %>%
          select(-ID)
        
        while (maxScore < threshold) {
          cons <- cons %>%
            select(-Score) %>%
            .[,-ncol(.)] %>%
            split(1:nrow(.)) %>%
            map(~{
              mutate(.,Score = prod(.,na.rm = T))  
            }) %>%
            bind_rows()
          
          maxScore <- max(cons$Score)
        }
        
        cons <- cons %>%
          mutate(ID = 1:nrow(.)) %>%
          filter(Score == max(Score)) %>%
          .[1,]
        
        consensusLevels <- names(cons)[1:(ncol(cons) - 2)]
        
        consensusClass <- classes %>%
          filter(ID == cons$ID) %>%
          select(MF:Adduct,consensusLevels) %>%
          mutate(Score = cons$Score)
        
        return(list(classes = classes,consensusScores = consensus,consensusClass = consensusClass))
      }) %>%
      map(~{
        .$consensusClass
      }) %>%
      bind_rows()
  } else {
    consensusClasses <- classifications %>%
      select(-CID,-InChIKey) %>%
      mutate(Score = 1)
  }
  return(consensusClasses) 
}

#' consensusClassification
#' @description Calculate a consensus classification for a given molecular formula and adducts.
#' @param MF molecular formula
#' @param adducts character vector of adducts
#' @param threshold consensus threshold
#' @examples
#' consensusClassification('C10H10O7')
#' @importFrom stringr str_detect
#' @importFrom tibble tibble
#' @export

consensusClassification <- function(MF, adducts = c('[M-H]1-'), threshold = 0.5){
  hits <- pubchemMatch(MF)
  
  if (is.null(hits)) {
    hits <- tibble()
    PIPs <- tibble()
    classifications <- tibble()
    con <-  tibble(MF = MF,Adduct = adducts,kingdom = 'No hits')
  } else {
    PIPs <- pips(hits,adducts)
    
    if (nrow(PIPs) == 0) {
      PIPs <- tibble()
      classifications <- tibble()
      con <- tibble(MF = MF,Adduct = adducts,kingdom = 'No hits')
    } else {
      classifications <- pipClassifications(PIPs)
      
      if (nrow(classifications) == 0) {
        classifications <- tibble()
        con <- tibble(MF = MF,Adduct = adducts,kingdom = 'Unclassified')
      } else {
        con <- classifications %>%
          consensusCls(threshold = threshold)   
      }
      
    }
  }
  
  consensus <- new('Consensus')
  consensus@hits <- hits
  consensus@PIPs <- PIPs
  consensus@classifications <- classifications
  consensus@consensus <- con
  
  return(consensus)
}

#' consensus
#' @rdname consensus
#' @description Consensus classifications for molecular formula assignments.
#' @param x S4 object of class Workflow
#' @param organism organism kegg ID.
#' @param threshold majority assignment threshold for consensus classifications
#' @importClassesFrom MFassign Assignment
#' @importFrom MFassign assignments
#' @importFrom methods new
#' @importFrom metaboWorkflows resultsAnnotation flags
#' @importFrom lubridate seconds_to_period
#' @importClassesFrom metaboWorkflows Workflow
#' @export

setMethod('consensus',signature = 'Workflow',
          function(x,organism = 'hsa', threshold = 0.5){
            
            if (!('MFassignment' %in% flags(x))) {
              stop('Molecular assignment needs to be run on Workflow object!')
            }
            
            a <- x %>%
              resultsAnnotation()
            
            z <- keggConsensus(a,organism = 'bdi')
            
            n <- z %>%
              filter(kingdom == 'No hits' | kingdom == 'Unclassified') %>%
              distinct()
            
            startTime <- proc.time()
            pc <- n %>%
              select(MF,Adduct) %>%
              split(.$MF) %>%
              map(~{
                consensusClassification(.$MF[1],.$Adduct)
              }) 
            endTime <- proc.time()
            
            elapsed <- {endTime - startTime} %>%
              .[3] %>%
              round(1) %>%
              seconds_to_period() %>%
              str_c('[',.,']')
            
            cat('\n',elapsed)
            
            return(pc)
          }
)