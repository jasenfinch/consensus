
setMethod('show','Consensus',
          function(object){
            cat('Consensus structural classifications\n\n')
            cat('MF:\t\t',mf(object),'\n')
            cat('Adducts:\t',nrow(adductRules(object)),'\n')
            cat('Organism:\t',organism(object),'\n')
            cat('Database:\t',database(object),'\n')
            cat('Threshold:\t',str_c(threshold(object),'%'),'\n')
            cat('Hits:\t\t',hits(object) %>% getAccessions() %>% nrow(),'\n')
            cat('Classifications:',nrow(classifications(object)),'\n')
            cat('PIPs:\t\t',nrow(PIPs(object)),'\n')
            cat('Consensuses:\t',nrow(consensusClassifications(object)))
          }
)