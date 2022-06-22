
setClass('Consensus',
         slots = list(
           MF = 'character',
           adduct_rules = 'tbl_df',
           organism = 'character',
           database = 'character',
           threshold = 'numeric',
           classifications = 'tbl_df',
           PIPs = 'tbl_df',
           consensus = 'tbl_df'
         ),
         contains = 'MetaboliteDatabase',
         prototype = list(
           classifications = tibble(),
           PIPs = tibble(),
           consensus = tibble()
         )
)

#' @importFrom cli tree
#' @importFrom stats na.omit

setGeneric('classificationTree',function(x){
  standardGeneric('classificationTree')
})

setMethod('classificationTree',signature = 'Consensus',
          function(x){
            d <- x %>%
              overallConsensus() 
            
            if (nrow(d) > 0){
              d <- d %>%
                select(-`Consensus (%)`) %>%
                gather(Level,Name) %>%
                mutate(Label = str_c(Level,': ',Name)) %>%
                na.omit() %>%
                select(Label)
              
              connections <- list()
              for (i in 1:nrow(d)) {
                if (i == nrow(d)) {
                  connections[[i]] <- c(character(0))
                } else {
                  connections[[i]] <- c(d$Label[(i + 1)])  
                }
              } 
              
              a <- data.frame(stringsAsFactors = FALSE,
                              id = d$Label,
                              connections = I(connections))
              
              tree(a)   
            } else {
              invisible()
            }
            
          }
)

setMethod('show','Consensus',
          function(object){
            cat('Consensus structural classifications\n\n')
            cat('MF:\t\t\t',mf(object),'\n')
            cat('Adducts:\t\t',nrow(adductRules(object)),'\n')
            cat('Organism:\t\t',organism(object),'\n')
            cat('Database:\t\t',database(object),'\n')
            cat('Threshold:\t\t',str_c(threshold(object),'%'),'\n')
            cat('Hits:\t\t\t',hits(object) %>% entries() %>% nrow(),'\n')
            cat('Classifications:\t',nrow(classifications(object)),'\n')
            cat('Average PIPs:\t\t',
                {
                  pips <- object %>%
                    PIPs()
                  
                  if (nrow(pips) > 0){
                    pips %>% 
                      group_by(Adduct) %>% 
                      summarise(N = n()) %>% 
                      .$N %>% 
                      mean() %>% 
                      round()
                  } else {
                    0
                  }
                },
                '\n')
            cat('Average Consensus:\t',str_c(
              {
                con <- object %>% 
                  consensusClassifications() 
                
                if (nrow(con) > 0){
                  con %>% 
                    .$`Consensus (%)` %>% 
                    mean() %>% 
                    round()
                } else {
                  0
                }
              },
              '%\n\n'))
            
            classification_tree <- classificationTree(object)
            
            if (!is.null(classification_tree)) print(classification_tree)
          }
)

#' Accessor methods for the `Consensus` and `Construction` S4 classes
#' @rdname access
#' @description  Accessor methods for the `Consensus` and `Construction` S4 classes.
#' @param x object of S4  class `Consensus` or `Construction`
#' @export

setGeneric('mf',function(x){
  standardGeneric('mf')
})

#' @rdname access

setMethod('mf',signature = 'Consensus',
          function(x){
            x@MF
          })

#' @rdname access
#' @export

setGeneric('adductRules',function(x){
  standardGeneric('adductRules')
})

#' @rdname access

setMethod('adductRules',signature = 'Consensus',
          function(x){
            x@adduct_rules
          })

#' @rdname access
#' @export

setGeneric('organism',function(x){
  standardGeneric('organism')
})

#' @rdname access

setMethod('organism',signature = 'Consensus',
          function(x){
            x@organism
          })

#' @rdname access
#' @export

setGeneric('database',function(x){
  standardGeneric('database')
})

#' @rdname access

setMethod('database',signature = 'Consensus',
          function(x){
            x@database
          })

#' @rdname access
#' @export

setGeneric('threshold',function(x){
  standardGeneric('threshold')
})

#' @rdname access

setMethod('threshold',signature = 'Consensus',
          function(x){
            x@threshold
          })

#' @rdname access
#' @export

setGeneric('hits',function(x){
  standardGeneric('hits')
})

#' @rdname access
#' @importFrom methods as

setMethod('hits',signature = 'Consensus',
          function(x){
            as(x,'MetaboliteDatabase')
          })

#' @rdname access
#' @export

setGeneric('PIPs',function(x){
  standardGeneric('PIPs')
})

#' @rdname access

setMethod('PIPs',signature = 'Consensus',
          function(x){
            x@PIPs
          })

#' @rdname access
#' @export

setGeneric('classifications',function(x){
  standardGeneric('classifications')
})

#' @rdname access

setMethod('classifications',signature = 'Consensus',
          function(x){
            x@classifications
          })

#' @rdname access
#' @export

setGeneric('consensusClassifications',function(x){
  standardGeneric('consensusClassifications')
})

#' @rdname access

setMethod('consensusClassifications',signature = 'Consensus',
          function(x){
            x@consensus
          })

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

setGeneric('consensus',function(x){
  standardGeneric('consensus')
})

setMethod('consensus',signature = 'Consensus',
          function(x){
            
            thresh <- threshold(x)
            p <- PIPs(x)
            
            if (nrow(p) > 0) {
              
              noPIPs <- adductRules(x) %>%
                select(Adduct = Name) %>%
                anti_join(p, by = "Adduct") %>%
                mutate(kingdom = 'No database hits',`Consensus (%)` = 100)
              
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
                                         kingdom = 'No database hits',
                                         `Consensus (%)` = 100)
            }
            
            x@consensus <- consensusClasses
            return(x) 
          }
)

setGeneric('overallConsensus',function(x){
  standardGeneric('overallConsensus')
})

setMethod('overallConsensus',signature = 'Consensus',
          function(x){
            con <- consensusClassifications(x) 
            
            if (nrow(con) > 0){
              con %>%
                select(-`Consensus (%)`) %>%
                conse(thresh = threshold(x))
            } else {
              con
            }
          })
