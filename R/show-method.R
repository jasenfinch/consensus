
setMethod('show','Consensus',
          function(object){
            cat('Consensus structural classifications\n\n')
            cat('MF:\t\t',mf(object),'\n')
            cat('Adducts:\t',nrow(adductRules(object)),'\n')
            cat('Organism:\t',organism(object),'\n')
            cat('Database:\t',database(object),'\n')
            cat('Threshold:\t',str_c(threshold(object),'%'),'\n')
            cat('Hits:',nrow(hits(object)),'\n')
            cat('PIPs:',nrow(PIPs(object)),'\n')
            cat('Classifications:',nrow(classifications(object)),'\n')
            cat('Consensuses:',nrow(consensusClassifications(object)))
          }
)