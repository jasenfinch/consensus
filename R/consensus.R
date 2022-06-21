#' Consensus
#' @description An S4 class to store consensus structral classification results for a molecular formula.
#' @slot MF molecular formula
#' @slot adduct_rules tibble containing adduct rules as returned by \code{mzAnnotation::adducts()}
#' @slot organism organism KEGG ID. NA if database is pubchem.
#' @slot database database, kegg or pubchem
#' @slot threshold \% majority for consensus classification
#' @slot classifications structural classifications for database hits
#' @slot PIPs putative ionisation product matches from database hits for given adduct rules
#' @slot consensus consensus structural classifications for adducts
#' @export

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
         contains = 'MetaboliteDatabase'
)

#' mf
#' @rdname mf
#' @description Get and the molecular formula of a Consensus object. 
#' @param x S4 object of class Consensus
#' @export

setGeneric('mf',function(x){
  standardGeneric('mf')
})

#' @rdname mf

setMethod('mf',signature = 'Consensus',
          function(x){
            x@MF
          })

#' adductRules
#' @rdname adductRules
#' @description Get the adduct rules of a Consensus object. 
#' @param x S4 object of class Consensus
#' @export

setGeneric('adductRules',function(x){
  standardGeneric('adductRules')
})

#' @rdname adductRules

setMethod('adductRules',signature = 'Consensus',
          function(x){
            x@adductRules
          })

#' organism
#' @rdname organism
#' @description Get the organism ID of a Consensus object. 
#' @param x S4 object of class Consensus
#' @export

setGeneric('organism',function(x){
  standardGeneric('organism')
})

#' @rdname organism

setMethod('organism',signature = 'Consensus',
          function(x){
            x@organism
          })

#' database
#' @rdname database
#' @description Get the database of a Consensus object. 
#' @param x S4 object of class Consensus
#' @export

setGeneric('database',function(x){
  standardGeneric('database')
})

#' @rdname database

setMethod('database',signature = 'Consensus',
          function(x){
            x@database
          })

#' threshold
#' @rdname threshold
#' @description Get the threshold of a Consensus object. 
#' @param x S4 object of class Consensus
#' @export

setGeneric('threshold',function(x){
  standardGeneric('threshold')
})

#' @rdname threshold

setMethod('threshold',signature = 'Consensus',
          function(x){
            x@threshold
          })

#' hits
#' @rdname hits
#' @description Get the database molecular formula matches of a Consensus object. 
#' @param x S4 object of class Consensus
#' @export

setGeneric('hits',function(x){
  standardGeneric('hits')
})

#' @rdname hits

setMethod('hits',signature = 'Consensus',
          function(x){
            x@hits
          })

#' PIPs
#' @rdname PIPs
#' @description Get the putative ionisation products of a Consensus object. 
#' @param x S4 object of class Consensus
#' @export

setGeneric('PIPs',function(x){
  standardGeneric('PIPs')
})

#' @rdname PIPs

setMethod('PIPs',signature = 'Consensus',
          function(x){
            x@PIPs
          })

#' classifications
#' @rdname classifications
#' @description Get the classifications of a Consensus object. 
#' @param x S4 object of class Consensus
#' @export

setGeneric('classifications',function(x){
  standardGeneric('classifications')
})

#' @rdname classifications

setMethod('classifications',signature = 'Consensus',
          function(x){
            x@classifications
          })

#' consensusClassifications
#' @rdname consensusClassifications
#' @description Get the consensus classifications of a Consensus object. 
#' @param x S4 object of class Consensus
#' @export

setGeneric('consensusClassifications',function(x){
  standardGeneric('consensusClassifications')
})

#' @rdname consensusClassifications

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

setGeneric('overallConsensus',function(x){
  standardGeneric('overallConsensus')
})

setMethod('overallConsensus',signature = 'Consensus',
          function(x){
            consensusClassifications(x) %>%
              select(-`Consensus (%)`) %>%
              conse(thresh = threshold(x))
          })
