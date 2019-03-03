
pipClassifications <- function(pips){
  classi <- pips %>%
    # filter(MolecularFormula == 'C4H6O5') %>%
    .$InChIKey %>%
    classify()
  
  classifications <- pips %>%
    left_join(classi, by = "InChIKey") %>%
    select(CID:Adduct,InChIKey,kingdom,superclass,class,subclass,`level 5`:names(.)[length(names(.))]) #%>%
  # filter(!is.na(kingdom))
  classifications$kingdom[is.na(classifications$kingdom)] <- 'Unclassified'
  
  return(classifications)
}