#' @importFrom tibble rowid_to_column
#' @importFrom dplyr everything group_by summarise right_join n
#' @importFrom tidyr gather
#' @importFrom tidyselect contains

consensus <- function(classifications,threshold = 0.5){
  
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
            mutate(.,Score = prod(.,na.rm = T) * 100)  
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
        
        return(consensusClass)
      }) %>%
      bind_rows()
  } else {
    consensusClasses <- classifications %>%
      select(-ACCESSION_ID,-InChI,-SMILES,-InChIKey) %>%
      mutate(Score = 1)
  }
  
  consensusClasses <- consensusClasses %>%
    rename(`Consensus (%)` = Score) %>%
    mutate(`Consensus (%)` = `Consensus (%)`)
  
  return(consensusClasses) 
}