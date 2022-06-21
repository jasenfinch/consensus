
globalVariables(c('Compound','Consensus (%)','Enzyme','InChI','SMILES','ID','Level','Label','INCHIKEY'))

#' construction
#' @description Build or add to and load a consensus classification library. 
#' @param MFs Molecular formulas and adducts to search. Should be a tibble containing two character columns named MF and Adduct.
#' @param library_path target file library_path for classification library for storing consensus classifications
#' @param db databases to search. Can be either kegg or pubchem.
#' @param organism KEGG organism ID. Ignored if kegg is not specified in db.
#' @param threshold percentage majority threshold for consensus classifications
#' @param adduct_rules_table data frame containing adduct formation rules. The defaults is `mzAnnotation::adduct_rules()`.
#' @param classyfireR_cache file library_path for a `classyfireR` cache. See the documentation of `classyfireR::get_classification` for more details. 
#' @examples 
#' \dontrun{
#' MFs <- tibble(MF = c(rep('C12H22O11',2),'C4H6O5'),
#'               Adduct = c('[M-H]1-','[M+Cl]1-','[M-H]1-'))
#' structural_classifications <- construction(MFs)
#' } 
#' @importFrom purrr walk
#' @export

construction <- function(MFs, 
                         library_path = tempdir(), 
                         db = c('kegg','pubchem'), 
                         organism = character(), 
                         threshold = 50,
                         adduct_rules_table = adduct_rules(),
                         classyfireR_cache = NULL){
  
  if (ncol(MFs) != 2 & !identical(names(MFs),c('MF','Adduct'))) {
    stop('Argument MFs should be a tibble containing two character columns named MF and Adduct')
  }
  
  if (length(organism) == 0) {
    org <- 'none'
  } else {
    org <- organism
  }
  
  db <- match.arg(db,
                  choices = c('kegg',
                              'pubchem'),
                  several.ok = TRUE) %>% 
    sort()
  
  mfs <- MFs$MF %>%
    map(~{
      tibble(MF = .,database = db)
    }) %>%
    bind_rows()
  
  library_path <- normalizePath(library_path)
  
  libraryPath <- paste0(library_path,'/','structural_classification_library')
  
  classificationLibrary <- loadLibrary(MFs,libraryPath)
  
  if (length(classificationLibrary) > 0){
    statuses <- classificationLibrary %>%
      map(status) %>%
      bind_rows()
    
    if (nrow(statuses) > 0) {
      statuses <- statuses %>%
        filter(database %in% db)
    }
    
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
  } else {
    toDo <- mfs
  }
  
  message(str_c(length(unique(toDo$MF))),' MFs to retrieve out of ',length(unique(mfs$MF)))
  
  toDo %>%
    split(.$MF) %>%
    walk(~{
      d <- .
      dbase <- d$database
      
      for (i in dbase){
        consense <- construct(d$MF[1],
                              db = i,
                              organism = organism,
                              threshold = threshold,
                              adduct_rules_table = adduct_rules_table,
                              classyfireR_cache = classyfireR_cache)
        
        saveConsensus(consense,library_path = libraryPath) 
        
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
  
  message('\nComplete!')
  
  classificationLibrary <- suppressMessages(loadLibrary(MFs,libraryPath))
  
  statuses <- classificationLibrary %>%
    map(status) %>%
    bind_rows() %>%
    rowid_to_column(var = 'ID') %>%
    filter(database %in% db,MF %in% {MFs$MF %>% 
        unique()})
  
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
  
  if (identical(db,c('kegg','pubchem'))) {
    statuses <- statuses %>%
      split(.$MF) %>%
      map(~{
        d <- .
        
        if (length(unique(d$database)) > 1) {
          if (d$status[d$database == 'kegg'] == 'Classified') {
            d <- d %>%
              filter(database == 'kegg')
          } else {
            if (d$status[d$database == 'kegg'] == 'No hits') {
              d <- d %>%
                filter(database == 'pubchem')
            } else {
              if ((d$status[d$database == 'kegg'] == 'Unclassified') & (d$status[d$database == 'pubchem'] != 'Classified')) {
                d <- d %>%
                  filter(database == 'kegg')
              } else {
                d <- d %>%
                  filter(database == 'pubchem')
              }  
            }  
            
          }
          
        }
        
        return(d)
      }) %>%
      bind_rows()
  }
  
  MFs %>%
    left_join(statuses, by = "MF") %>%
    split(1:nrow(.)) %>%
    map(~{
      d <- .
      con <- classificationLibrary[[d$ID[1]]] %>%
        consensusClassifications() %>%
        filter(Adduct == d$Adduct[1])
      d %>%
        left_join(con, by = "Adduct") %>%
        select(-ID,-status)
    }) %>%
    bind_rows() %>%
    select(`Consensus (%)`,everything()) %>%
    select(MF:last_col(),`Consensus (%)`)
}

toSearch <- function(mfs_status,db){
  to <- mfs_status %>%
    split(.$MF) %>%
    map(~{
      if (!('Classified' %in% .x$status)) {
        return(.x)
      } else {
        tibble(MF = character(),
               database = character(),
               organism = character(),
               status = character())
      }
    }) %>%
    bind_rows() %>%
    filter(status != 'No hits' | is.na(status))
  
  if (nrow(to) > 0) {
    to <- to %>%
      split(.$MF) %>%
      map(~{
        if (('kegg' %in% .x$database)){
          .x <- .x %>%
            filter(database != 'pubchem')
        }
        return(.x)
      }) %>%
      bind_rows()
  }
  to %>%
    filter(status != 'Unclassified' | is.na(status)) %>%
    select(MF,database) %>%
    distinct()
}
