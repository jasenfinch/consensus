
pipClassifications <- function(PIPs,nCores = detectCores() * 0.75){
 
  inchis <- PIPs %>%
    .$InChIKey %>%
    unique() 
  
  message('Retreiving classifications...')
  
  slaves <- length(inchis) / 100
  slaves <-  ceiling(slaves)
  
  if (slaves > nCores) {
    slaves <- nCores
  }
  
  classi <- inchis %>% 
    classify(slaves)
  
  classifications <- PIPs %>%
    left_join(classi, by = "InChIKey") %>%
    select(CID,MolecularFormula,Adduct,InChIKey,kingdom,superclass,class,subclass,`level 5`:names(.)[length(names(.))]) #%>%
  # filter(!is.na(kingdom))
  
  classifications$kingdom[is.na(classifications$kingdom)] <- 'Unclassified'
  
  message(str_c(nrow(unique(classifications$InChIKey) %>% na.omit()),' classifications returned'))
  
  return(classifications)
}