#' @export

pipClassifications <- function(pips,nCores = detectCores(), clusterType = 'PSOCK'){
  classi <- pips %>%
    # filter(MolecularFormula == 'C4H6O5') %>%
    .$InChIKey %>%
    unique() 
  
  message(length(classi),' InChIKeys to classify')
  
  classi <- classi %>% 
    classify(nCores,clusterType)
  
  classifications <- pips %>%
    left_join(classi, by = "InChIKey") %>%
    select(CID:Adduct,InChIKey,kingdom,superclass,class,subclass,`level 5`:names(.)[length(names(.))]) #%>%
  # filter(!is.na(kingdom))
  classifications$kingdom[is.na(classifications$kingdom)] <- 'Unclassified'
  
  return(classifications)
}