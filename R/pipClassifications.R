
setMethod('pipClassifications',signature = 'Consensus',
          function(ips,nCores = detectCores(), clusterType = 'PSOCK'){
            pips <- ips@PIPs
            classi <- pips %>%
              .$InChIKey %>%
              unique() 
            
            message(length(classi),' InChIKeys to classify')
            
            classi <- classi %>% 
              classify(nCores,clusterType)
            
            classifications <- pips %>%
              left_join(classi, by = "InChIKey") %>%
              select(CID:Adduct,InChIKey,kingdom,superclass,class,subclass,`level 5`:names(.)[length(names(.))])
            classifications$kingdom[is.na(classifications$kingdom)] <- 'Unclassified'
            
            return(classifications)
          }
)