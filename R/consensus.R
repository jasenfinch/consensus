
consensus <- function(classifications){
  consensusClasses <- classifications %>%
   split(str_c(.$MolecularFormula,.$Adduct,sep = ' ')) %>%
    map(~{
      cl <- .
      levels <- names(cl)[5:length(names(cl))]
      
      classes <- cl %>%
        select(MolecularFormula,Adduct,kingdom:names(.)[length(names(.))]) %>%
        distinct() %>%
        rowid_to_column(var = 'ID')
      
      freq <- cl %>%
        left_join(classes, by = c("MolecularFormula", "Adduct", "kingdom", "superclass", "class", "subclass", "level 5", "level 6", "level 7", "level 8", "level 9")) %>%
        group_by(ID) %>%
        summarise(N = n()) %>%
        right_join(classes, by = "ID")
      
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
                group_by(MolecularFormula,Adduct,Class) %>%
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
            left_join(votes, by = c("MolecularFormula", "Adduct", "Level", "Class")) %>%
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
      
      while (maxScore < 0.5) {
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
        select(MolecularFormula:Adduct,consensusLevels) %>%
        mutate(Score = cons$Score)
      
      return(list(classes = classes,consensusScores = consensus,consensusClass = consensusClass))
    }) %>%
    map(~{
      .$consensusClass
    }) %>%
    bind_rows()
}