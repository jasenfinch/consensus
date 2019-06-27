#' @importFrom dplyr left_join
pipClassifications <- function(PIPs){
 
  inchis <- PIPs %>%
    .$InChIKey %>%
    unique() 
  
  message('Retreiving classifications...')
  
  classi <- inchis %>% 
    classify()
  
  classifications <- PIPs %>%
    left_join(classi, by = "InChIKey") %>%
    select(CID,MolecularFormula,Adduct,InChIKey,kingdom,superclass,class,subclass,`level 5`:names(.)[length(names(.))]) %>%
    filter(!is.na(kingdom))
  
  # classifications$kingdom[is.na(classifications$kingdom)] <- 'Unclassified'
  
  message(str_c(length(unique(classifications$InChIKey)),' classifications returned'))
  
  # classifications <- classifications %>%
    # filter(kingdom == 'Unclassified')
  
  return(classifications)
}