#' @importClassesFrom cheminf MetaboliteDatabase

setClass('Consensus',
         slots = list(
           MF = 'character',
           adduct_rules = 'tbl_df',
           organism = 'character',
           database = 'character',
           classifications = 'tbl_df',
           PIPs = 'tbl_df',
           consensus = 'tbl_df'
         ),
         contains = 'MetaboliteDatabase',
         prototype = list(
           classifications = tibble(
             ID = character(),
             INCHIKEY = character(),
             kingdom = character()
           ),
           PIPs = tibble(),
           consensus = tibble()
         )
)

#' @importFrom cli tree
#' @importFrom stats na.omit

setGeneric('classificationTree',function(x,adduct,threshold = 66){
  standardGeneric('classificationTree')
})

#' @importFrom dplyr any_of

setMethod('classificationTree',signature = 'Consensus',
          function(x,adduct,threshold = 66){
            
            if (nrow(classifications(x)) > 0){
              d <- x %>%
                consensus(adduct,threshold) %>%
                select(
                  any_of(c('kingdom',
                           'superclass',
                           'class',
                           'subclass')),
                  any_of(paste('level',5:10))) %>%
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

#' @importFrom dplyr summarise

setMethod('show','Consensus',
          function(object){
            cat('Consensus structural classifications\n\n')
            cat('MF:\t\t\t',mf(object),'\n')
            cat('Adducts:\t\t',nrow(adductRules(object)),'\n')
            cat('Organism:\t\t',organism(object),'\n')
            cat('Database:\t\t',database(object),'\n')
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
                '\n\n')
            
            classification_tree <- classificationTree(object,'M')
            
            if (!is.null(classification_tree)) print(classification_tree)
          }
)

#' Accessor methods for the `Consensus` and `Construction` S4 classes
#' @rdname access
#' @description  Accessor methods for the `Consensus` and `Construction` S4 classes.
#' @param x object of S4  class `Consensus` or `Construction`
#' @param adduct the ionisation adduct for which a consensus should be calculated 
#' @param threshold the percentage majority threshold for consensus classification
#' @details 
#' * `mf` - Return the searched molecular formula
#' * `adductRules` - Return a tibble of adduct formation rules.
#' * `organism` - Return the KEGG organism ID.
#' * `database` - Return the searched database.
#' * `threshold` - Return the percentage consensus threshold for structural classification selection.
#' * `hits` - Return a `MetaboliteDatabase` ionisation database of matched database hits.
#' * `PIPs` - Return the putative ionisation products of database hits
#' * `classifications` -Return the structural chemical classifications of database hits.
#' * `consensusClassifications` - Return the consensus classification or classifications.
#' * `summariseClassifications` - Return a tibble of summarised consensus structural classifications.
#' @return 
#' A character, a numeric, a tibble or an object of S4 class `MetaboliteDatabase`, depending on the method used.
#' @examples 
#' consensus <- construct(
#'   'C4H6O5',
#'   organism = 'hsa')
#' 
#' ## Return the molecular formula
#' mf(consensus)
#' 
#' ## Return the adduct formation rules
#' adductRules(consensus)
#' 
#' ## Return the KEGG organism ID
#' organism(consensus)
#' 
#' ## Return the searched database
#' database(consensus)
#' 
#' ## Return the `MetaboliteDatabase` ionisation database of searched database  hits
#' hits(consensus)
#' 
#' ## Return the putative ionisation products
#' PIPs(consensus)
#' 
#' ## Return the structural classifications
#' classifications(consensus)
#' 
#' ## Return the consensus structural classification
#' consensusClassifications(consensus)
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

#' @rdname access
#' @export

setGeneric('consensus',function(x,adduct,threshold = 66){
  standardGeneric('consensus')
})

#' @rdname access

setMethod('consensus',signature = 'Consensus',
          function(x,adduct,threshold = 66){
            
            selected_adduct <- adduct
            
            classifications <- x %>% 
              consensusClassifications() %>% 
              filter(
                adduct == selected_adduct
              )
            
            above_threshold <- classifications %>% 
              filter(
                `consensus (%)` >= threshold
              ) 
            
            if (nrow(above_threshold) > 0) {
              above_threshold %>% 
                group_by(level) %>% 
                slice(1) %>% 
                ungroup() %>% 
                select(-n) %>% 
                mutate(
                  `consensus (%)` = min(`consensus (%)`)
                ) %>% 
                spread(
                  level,
                  class
                )
            } else {
              tibble(
                adduct = adduct,
                total = nrow(classifications),
                `consensus (%)` = 100,
                kingdom = 'No consensus'
              )
            }
          })

#' @importFrom dplyr group_by n count join_by
#' @importFrom tidyr gather
#' @importFrom tidyselect contains

setGeneric('calcConsensus',function(x){
  standardGeneric('calcConsensus')
})

setMethod('calcConsensus',signature = 'Consensus',
          function(x){
            
            x@consensus <- x %>% 
              adductRules() %>% 
              select(adduct = Name) %>% 
              left_join(
                x %>% 
                  PIPs(),
                join_by(adduct == Adduct)) %>% 
              left_join(
                x %>% 
                  classifications() %>% 
                  mutate(
                    ID = as.numeric(ID)
                  ),
                join_by(ID)
              ) %>% 
              mutate(
                kingdom = kingdom %>% 
                  replace(
                    is.na(kingdom),
                    'No database hits'
                  )
              ) %>% 
              group_by(
                adduct
              ) %>% 
              mutate(
                total = n()
              ) %>% 
              gather(
                level,
                class,
                -(ID:INCHIKEY),
                -adduct,
                -total) %>% 
              filter(
                !is.na(class)
              ) %>% 
              group_by(
                adduct,
                level,
                class,
                total
              ) %>% 
              count() %>% 
              mutate(
                `consensus (%)` = n / total * 100,
                level = factor(
                  level,
                  levels = c(
                    'kingdom',
                    'superclass',
                    'class',
                    'subclass',
                    paste('level',5:10)
                  ))
              ) %>% 
              ungroup()
            
            return(x) 
          }
)
