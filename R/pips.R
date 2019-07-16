#' @importFrom dplyr left_join

pipClassifications <- function(PIPs){
 
  inchis <- PIPs %>%
    .$InChIKey %>%
    unique() 
  
  classi <- inchis %>% 
    classify()
  
  if (nrow(classi) > 0) {
    classifications <- PIPs %>%
      select(-(SMILE:ACCESSION_ID),CID,MF,Adduct,InChIKey) %>%
      left_join(classi, by = "InChIKey") %>%
      select(CID,MF,Adduct,InChIKey,kingdom,everything()) %>%
      filter(!is.na(kingdom))
    
    # classifications$kingdom[is.na(classifications$kingdom)] <- 'Unclassified'
    
    # classifications <- classifications %>%
    # filter(kingdom == 'Unclassified')
    
    return(classifications)   
  } else {
    return(classi)
  }
 
}
