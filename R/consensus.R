conse <- function(cl,thresh){
  thresh <- thresh / 100
  
  levels <- names(cl)[which(names(cl) == 'kingdom'):length(names(cl))]
  
  classes <- cl %>%
    select(-Adduct) %>%
    distinct() %>%
    rowid_to_column(var = 'id')
  
  suppressMessages(freq <- cl %>%
                     left_join(classes) %>%
                     group_by(id) %>%
                     summarise(N = n()) %>%
                     right_join(classes, by = "id"))
  
  votes <- levels %>%
    map(~{
      lev <- .
      freq %>%
        rename('Class' = !!lev) %>%
        select(id,Class,N) %>%
        split(.$Class) %>%
        map(~{
          d <- .
          d %>%
            group_by(Class) %>%
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
        gather('Level','Class',-id) %>%
        left_join(votes, by = c("Level", "Class")) %>%
        select(-Class) %>%
        spread(Level,N)
    }) %>%
    bind_rows() 
  
  clLevels <- clLevels[clLevels %in% names(votesTable)]
  
  votesTable <- votesTable %>%
    select(id,clLevels,contains('level'))
  
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
    mutate(id = 1:nrow(.))
  
  consensus <- proportions %>%
    select(-id) %>% 
    split(1:nrow(.)) %>%
    map(~{
      mutate(.,Score = prod(.,na.rm = T))  
    }) %>%
    bind_rows() %>%
    mutate(id = 1:nrow(.))
  
  maxScore <- max(consensus$Score)
  
  cons <- consensus %>%
    select(-id)
  
  while (maxScore < thresh) {
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
    mutate(id = 1:nrow(.)) %>%
    filter(Score == max(Score)) %>%
    .[1,]
  
  consensusLevels <- names(cons)[1:(ncol(cons) - 2)]
  
  consensusClass <- classes %>%
    filter(id == cons$id) %>%
    select(consensusLevels) %>%
    mutate(`Consensus (%)` = cons$Score * 100)
  
  return(consensusClass)
}

#' @importFrom tibble rowid_to_column
#' @importFrom dplyr everything group_by summarise right_join n anti_join full_join
#' @importFrom tidyr gather
#' @importFrom tidyselect contains

setMethod('consensus',signature = 'Consensus',
          function(x){
            
            thresh <- threshold(x)
            p <- PIPs(x)
            
            if (nrow(p) > 0) {
              
              noPIPs <- adductRules(x) %>%
                select(Adduct = Name) %>%
                anti_join(p, by = "Adduct") %>%
                mutate(kingdom = 'No hits',`Consensus (%)` = 100)
              
              classi <- classifications(x) %>%
                filter(kingdom != 'Unclassified') %>%
                filter(ID %in% p$ID) %>%
                left_join(p,by = c('ID')) %>%
                select(Adduct,everything())
              
              noClassi <- p %>%
                anti_join(classi, by = c("Adduct")) %>%
                select(-ID) %>%
                mutate(kingdom = 'Unclassified',`Consensus (%)` = 100) 
              
              if (nrow(classi) > 1) {
                consensusClasses <- classi %>%
                  split(.$Adduct) %>%
                  map(conse,thresh = thresh) %>%
                  bind_rows(.id = 'Adduct') %>%
                  full_join(noPIPs, by = c("Adduct", "kingdom", "Consensus (%)")) %>%
                  full_join(noClassi, by = c("Adduct", "kingdom", "Consensus (%)")) %>%
                  distinct()
              } else {
                consensusClasses <- classi %>%
                  select(-ID,-INCHIKEY) %>%
                  mutate(`Consensus (%)` = 100)  %>%
                  full_join(noPIPs, by = c("Adduct", "kingdom", "Consensus (%)")) %>%
                  full_join(noClassi, by = c("Adduct", "kingdom", "Consensus (%)")) %>%
                  distinct()
              }  
            } else {
              consensusClasses <- tibble(Adduct = adductRules(x)$Name,
                                        kingdom = 'No hits',
                                        `Consensus (%)` = 100)
            }
            
            x@consensus <- consensusClasses
            return(x) 
          }
)

setMethod('overallConsensus',signature = 'Consensus',
          function(x){
            consensusClassifications(x) %>%
              select(-`Consensus (%)`) %>%
              conse(thresh = threshold(x))
          })
