#' @importFrom cli tree

setMethod('classificationTree',signature = 'Consensus',
          function(x){
            d <- x %>%
              overallConsensus() %>%
              select(-`Consensus (%)`) %>%
              gather(Level,Name) %>%
              mutate(Label = str_c(Level,': ',Name)) %>%
              select(Label)
            
            connections <- list()
            for (i in 1:nrow(d)){
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
            cat('Hits:\t\t\t',hits(object) %>% getAccessions() %>% nrow(),'\n')
            cat('Classifications:\t',nrow(classifications(object)),'\n')
            cat('Average PIPs:\t\t',
                object %>%
                  PIPs() %>% 
                  group_by(Adduct) %>% 
                  summarise(N = n()) %>% 
                  .$N %>% 
                  mean() %>% 
                  round(),'\n')
            cat('Average Consensus:\t',str_c(object %>% 
                  consensusClassifications() %>% 
                  .$`Consensus (%)` %>% 
                  mean() %>% 
                  round(),'%\n\n'))
            print(classificationTree(object))
          }
)