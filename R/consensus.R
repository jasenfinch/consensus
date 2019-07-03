
globalVariables(c('.','kingdom','CID','MF','Adduct','InChIKey','superclass','subclass','level 5','MF',
                  'Charge','CanonicalSMILES','CovalentUnitCount'))

#' @importFrom tibble rowid_to_column
#' @importFrom dplyr everything group_by summarise right_join
#' @importFrom tidyr gather

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
          bind_rows() %>%
          select(ID:Adduct,kingdom,superclass,class,subclass,`level 5`:names(cl)[length(names(cl))])
        
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
    return(tibble(MF = MF,Adduct = adducts,kingdom = 'No hits'))
  }
  
  PIPs <- pips(hits,adducts)
  
  if (nrow(PIPs) == 0) {
    return(tibble(MF = MF,Adduct = adducts,kingdom = 'No hits'))
  }
  
  classifications <- pipClassifications(PIPs)
  
  if (nrow(classifications) == 0) {
    return(tibble(MF = MF,Adduct = adducts,kingdom = 'Unclassified'))
  }
  
  con <- classifications %>%
    consensusCls(threshold = threshold)
  
  consensus <- new('Consensus')
  consensus@hits <- hits
  consensus@PIPs <- PIPs
  consensus@classifications <- classifications
  consensus@consensus <- con
  
  return(consensus)
}

#' @importClassesFrom MFassign Assignment
#' @importFrom MFassign assignments
#' @importFrom methods new

setMethod('consensus',signature = 'Assignment',
          function(x,organism = 'hsa',threshold = 0.5){
            
            ips <- x %>%
              assignments() %>%
              select(MF,Adduct) %>%
              distinct()
            
            con <- new('Consensus')
            con@IPs <- ips
            
            message(str_c('\nIdentifying consensus classifications for ',nrow(ips),' ionisation products consisting of ',length(unique(ips$MF)),' molecular formulas'))
            
            con <- pubchemPIPs(con)
            classifications <- pipClassifications(pips)
            
            consensusClasses <- classifications %>%
              filter(kingdom != 'Unclassified') %>%
              consensusCls(threshold = 1/3) %>%
              select('MF','Adduct','Score','kingdom','superclass','class','subclass',names(.)[str_detect(names(.),'level')])
            
            return(consensusClasses)
          }
)