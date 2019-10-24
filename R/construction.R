
globalVariables(c('Compound','Consensus (%)','Enzyme','InChI','SMILES','ID','Level','Label','INCHIKEY'))

#' construction
#' @description Build or add to and load a consensus classification library. 
#' @param MFs Molecular formulas and adducts to search. Should be a tibble containing two character columns named MF and Adduct.
#' @param path target file path for classification library 
#' @param db databases to search. Can be either kegg or pubchem.
#' @param organism KEGG organism ID. Ignored if kegg is not specified in db.
#' @param adductRules adduct rules table as returned by mzAnnotation::adducts()
#' @param threshold \% majority threshold for consensus classifications
#' @examples 
#' \dontrun{
#' MFs <- tibble(MF = c(rep('C12H22O11',2),'C4H6O5'),Adduct = c('[M-H]1-','[M+Cl]1-'))
#' structural_classifications <- construction()
#' } 
#' @importFrom purrr walk
#' @export

construction <- function(MFs, path = '.', db = c('kegg','pubchem'), organism = character(), threshold = 50){
  
  if (ncol(MFs) != 2 & !identical(names(MFs),c('MF','Adduct'))) {
    stop('Argument MFs should be a tibble containing two character columns named MF and Adduct')
  }
  
  if (F %in% ((MFs$Adduct %>% unique()) %in% mzAnnotation::adducts()$Name)) {
    stop('Adducts should be one of those available in mzAnnotation::adducts()$Name')
  }
  
  if (length(organism) == 0) {
    org <- 'none'
  } else {
    org <- organism
  }
  
  mfs <- MFs$MF %>%
    map(~{
      tibble(MF = .,database = db)
    }) %>%
    bind_rows()
  
  libraryPath <- str_c(path,'classification_library',sep = '/')
  
  if (isFALSE(checkLibrary(path))) {
    dir.create(libraryPath)  
    message(str_c('\nCreated structural classifcation library at ',libraryPath))
    toDo <- mfs
  } else {
    message(str_c('\nUsing structural classifcation library at ',libraryPath))
    
    classificationLibrary <- loadLibrary(path)
    
    statuses <- classificationLibrary %>%
      map(status) %>%
      bind_rows() %>%
      filter(database %in% db)
    
    if ('kegg' %in% db) {
      statuses <- statuses %>%
        split(.$database) %>%
        map(~{
          d <- .
          if (d$database[1] == 'kegg') {
            d <- d %>% 
              filter(organism == org)  
          }
          return(d)
        }) %>%
        bind_rows()
    }
    
    mfs_status <- mfs %>%
      left_join(statuses, by = c("MF", "database"))
    
    toDo <- mfs_status %>%
      toSearch(db)
  }
    
    message(str_c(length(unique(toDo$MF))),' MFs to retrieve out of ',length(unique(mfs$MF)))
    
    toDo %>%
      split(.$MF) %>%
      walk(~{
        d <- .
        dbase <- d$database
        
        for (i in dbase){
          consense <- construct(d$MF[1],db = i,organism = organism,threshold = threshold)
          
          saveConsensus(consense,path = libraryPath) 
          
          if (database(consense) == 'kegg') {
            kingdoms <- consense %>%
              consensusClassifications() %>%
              .$kingdom %>%
              unique() %>%
              sort()
            
            if (!identical(kingdoms,'No hits') & !identical(kingdoms,'Unclassified') & !identical(kingdoms,c('No hits','Unclassified'))) {
              break()
            }
          }
        }
      })
    
    classificationLibrary <- loadLibrary(path)
    
    statuses <- classificationLibrary %>%
      map(status) %>%
      bind_rows() %>%
      filter(database %in% db)
    
    if ('kegg' %in% db) {
      statuses <- statuses %>%
        split(.$database) %>%
        map(~{
          d <- .
          if (d$database[1] == 'kegg') {
            d <- d %>% 
              filter(organism == org)  
          }
          return(d)
        }) %>%
        bind_rows()
    }
    
    mfs_status <- mfs %>%
      left_join(statuses, by = c("MF", "database"))
  
  message('\nComplete!')
  
}

toSearch <- function(mfs,db){
  mfs %>%
    filter(is.na(status)) %>%
    split(.$MF) %>%
    map(~{
      d <- .
      if (('kegg' %in% db) & !('kegg' %in% d$database)){
        d <- d %>%
          filter(database != 'pubchem')
      }
      return(d)
    }) %>%
    bind_rows() %>%
    select(MF,database)
}
