
globalVariables(c('Compound','Consensus (%)','Enzyme','InChI','SMILES','ID','Level','Label',
                  'INCHIKEY','.','Adduct','CanonicalSMILES','Charge','CovalentUnitCount',
                  'Feature','INCHI','IUPACName','InChIKey','Isotope','MF','MolecularFormula',
                  'Name','Score','id','kingdom','level','value'))

setClass('Construction',
         slots = list(
           classifications = 'tbl_df'
         ),
         contains = 'Assignment',
         prototype = list(
           classifications = tibble()
         ))

setMethod('show',signature = 'Construction',
          function(object){
            cat('Consensus structural classifications','\n')
            cat(paste0('Assignments: ',nrow(assignments(object))),'\n')
            cat(paste0('Classifications: ',nrow(object@classifications)))
          })

#' @rdname access
#' @importFrom assignments featureData

setMethod('classifications',signature = 'Construction',
          function(x) {
            x %>% 
              featureData() %>%
              {tibble(Feature = colnames(.))} %>% 
              left_join(assignments(x),
                        by = 'Feature') %>% 
              select(Feature,Name,MF,Isotope,Adduct) %>% 
              left_join(x@classifications,
                        by = c('MF','Adduct')) %>% 
              mutate(kingdom = kingdom %>% 
                       replace(is.na(MF),
                               'Unknown'))
            })

#' @rdname access
#' @export

setGeneric('summariseClassifications',function(x)
  standardGeneric('summariseClassifications'))

#' @rdname access
#' @importFrom tidyr drop_na

setMethod('summariseClassifications',signature = 'Construction',
          function(x){
            x %>% 
              classifications() %>% 
              select(kingdom:last_col(1)) %>% 
              gather(level,value) %>% 
              drop_na() %>% 
              group_by(level,value) %>% 
              summarise(count = n(),
                        .groups = 'drop')
          })

#' Consensus structural classifications
#' @rdname construction
#' @description Build or add to and load a consensus classification library. 
#' @param x Molecular formulas and adducts to search. Should either be a tibble containing two character columns named MF and Adduct or and S4 object of class `Assignment`.
#' @param library_path target file library_path for classification library for storing consensus classifications
#' @param db databases to search. Can be either `kegg` and/or `pubchem`.
#' @param organism KEGG organism ID. Ignored if kegg is not specified in db.
#' @param threshold percentage majority threshold for consensus classifications
#' @param adduct_rules_table data frame containing adduct formation rules. The defaults is `mzAnnotation::adduct_rules()`.
#' @param classyfireR_cache file library_path for a `classyfireR` cache. See the documentation of `classyfireR::get_classification` for more details. 
#' @return If argument `x` is a tibble, then a tibble is returned containing the consensus structural classifications. If argument `x` is an object of S4 class `Assignment`, and object of S4 class `Construction` is returned.
#' @examples 
#' \dontrun{
#' x <- tibble::tibble(MF = c(rep('C12H22O11',2),'C4H6O5'),
#'               Adduct = c('[M-H]1-','[M+Cl]1-','[M-H]1-'))
#' structural_classifications <- construction(x)
#' } 
#' @importFrom purrr walk2
#' @importFrom tidyr expand_grid
#' @export

setGeneric('construction',function(x, 
                                   library_path = tempdir(), 
                                   db = 'kegg', 
                                   organism = character(), 
                                   threshold = 50,
                                   adduct_rules_table = adduct_rules(),
                                   classyfireR_cache = NULL)
  standardGeneric('construction')
)

#' @rdname construction

setMethod('construction',signature = 'tbl_df',
          function(x, 
                   library_path = tempdir(), 
                   db = 'kegg', 
                   organism = character(), 
                   threshold = 50,
                   adduct_rules_table = adduct_rules(),
                   classyfireR_cache = NULL){
            
            if (ncol(x) != 2 & !identical(names(x),c('MF','Adduct'))) {
              stop('Argument x should be a tibble containing two character columns named MF and Adduct')
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
            
            mfs <- expand_grid(
              MF = x$MF,
              database = db
            )
            
            library_path <- normalizePath(library_path)
            
            libraryPath <- paste0(library_path,'/','construction_library')
            
            classificationLibrary <- loadLibrary(x,libraryPath)
            
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
                  filter(database != 'kegg' | 
                           (database == 'kegg' & organism == org))
              }
              
              mfs_status <- mfs %>%
                left_join(statuses, 
                          by = c("MF", "database"))
              
              toDo <- mfs_status %>%
                toSearch(db)
            } else {
              toDo <- mfs
            }
            
            message(str_c(length(unique(toDo$MF))),' MFs to retrieve out of ',length(unique(mfs$MF)))
            
            search_mfs <- toDo %>%
              rowid_to_column(var = 'idx') %>% 
              split(.$MF) 
            
            n_mfs <- length(search_mfs)
            
            search_mfs %>%
              walk2(seq_along(.),
                    ~{
                
                message(' ')
                message(paste0(.y,'. (',round(.y / n_mfs * 100),'%) '),appendLF = FALSE)
                
                for (i in .x$database){
                  consense <- construct(.x$MF[1],
                                        db = i,
                                        organism = organism,
                                        threshold = threshold,
                                        adduct_rules_table = adduct_rules_table,
                                        classyfireR_cache = classyfireR_cache)
                  
                  saveConsensus(consense,path = libraryPath) 
                  
                  if (database(consense) == 'kegg') {
                    kingdoms <- consense %>%
                      consensusClassifications() %>%
                      .$kingdom %>%
                      unique() %>%
                      sort()
                    
                    if (!identical(kingdoms,'No database hits') & 
                        !identical(kingdoms,'Unclassified') & 
                        !identical(kingdoms,c('No database hits','Unclassified'))) {
                      break()
                    }
                  }
                }
              })
            
            message('\nComplete!')
            
            classificationLibrary <- suppressMessages(loadLibrary(x,libraryPath))
            
            statuses <- classificationLibrary %>%
              map(status) %>%
              bind_rows() %>%
              rowid_to_column(var = 'ID') %>%
              filter(database %in% db,MF %in% {x$MF %>% 
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
                      if (d$status[d$database == 'kegg'] == 'No database hits') {
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
            
            x %>%
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
)

#' @importFrom assignments assignments
#' @rdname construction

setMethod('construction',signature = 'Assignment',
          function(x, 
                   library_path = tempdir(), 
                   db = 'kegg', 
                   organism = character(), 
                   threshold = 50,
                   classyfireR_cache = NULL){
            
            adduct_rules_table <- assignments::adductRules(x)
            
            mfs <- x %>% 
              assignments() %>% 
              select(MF,Adduct) %>% 
              distinct()
            
            structural_classifications <- construction(
              mfs,
              library_path = library_path,
              db = db,
              organism = organism,
              threshold = threshold,
              adduct_rules_table = adduct_rules_table,
              classyfireR_cache = classyfireR_cache
            )
            
            new('Construction',
                x,
                classifications = structural_classifications)
          }
)

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
    filter(status != 'No database hits' | is.na(status))
  
  to %>%
    filter(status != 'Unclassified' | is.na(status)) %>%
    select(MF,database) %>%
    distinct()
}

