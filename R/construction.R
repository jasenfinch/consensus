
globalVariables(c('Compound','consensus (%)','Enzyme','InChI','SMILES','ID','Level','Label',
                  'INCHIKEY','.','Adduct','CanonicalSMILES','Charge','CovalentUnitCount',
                  'Feature','INCHI','IUPACName','InChIKey','Isotope','MF','MolecularFormula',
                  'Name','Score','id','kingdom','level','value','file_name','adduct','total'))

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
                        by = c('MF' = 'mf','Adduct' = 'adduct')) %>% 
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

#' Consensus structural classifications for putative ionisation products
#' @rdname construction
#' @description Perform consensus structural classification for molecular formulas assigned to *m/z* features from electrospray ionisation mass spectrometry approaches. 
#' @param x The molecular formulas and adducts to search. This should either be a tibble containing two character columns named `MF` and `Adduct` or and S4 object of class `Assignment`.
#' @param library_path the target file path for the classification library in which to store consensus classification data
#' @param db the databases to search. This can either be `kegg` and/or `pubchem`.
#' @param organism the KEGG organism ID. This is Ignored if argument `db` is set to `pubchem` 
#' @param threshold the percentage majority threshold for consensus classification
#' @param adduct_rules_table a data frame containing the adduct formation rules. The defaults is `mzAnnotation::adduct_rules()`.
#' @param classyfireR_cache the file path for a `classyfireR` cache. See the documentation of `classyfireR::get_classification` for more details. 
#' @return If argument `x` is a tibble, then a tibble is returned containing the consensus structural classifications. If argument `x` is an object of S4 class `Assignment`, and object of S4 class `Construction` is returned.
#' @examples 
#' x <- tibble::tibble(MF = c('C12H22O11','C4H6O5'),
#'                     Adduct = c('[M+Cl]1-','[M-H]1-'))
#' structural_classifications <- construction(x)
#' 
#' structural_classifications
#' @importFrom purrr map_dfr
#' @importFrom dplyr cross_join group_split slice relocate arrange
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
            
            db <- match.arg(db,
                            choices = c('kegg',
                                        'pubchem'),
                            several.ok = TRUE) %>% 
              sort()
            
            if (length(organism) == 0){
              organism <- 'none'
            }
            
            library_path <- normalizePath(library_path) %>% 
              paste0(.,'/','construction_library')
            
            items <- cross_join(
              x,
              tibble(database = db)
            ) %>% 
              mutate(
                organism = organism
              )
            
            search_mfs <- items %>% 
              status(library_path = library_path) %>% 
              searchNecessary(db = db)
            
            message(str_c(length(unique(search_mfs$MF))),' MFs to retrieve out of ',length(unique(x$MF)))
            
            while(nrow(search_mfs) > 0){
              message()
              
              construct(
                MF = search_mfs$MF[1],
                db = search_mfs$database[1],
                organism = if (search_mfs$organism[1] == 'none') {character()} else {search_mfs$organism[1]},
                adduct_rules_table = adduct_rules_table,
                classyfireR_cache = classyfireR_cache
              ) %>% 
                saveConsensus(path = library_path) 
              
              search_mfs <- items %>% 
                status(library_path = library_path) %>% 
                searchNecessary(db = db)
            }
            
            message('\nComplete!')
            
            statuses <- items %>%
              status(library_path = library_path) %>% 
              filter(
                exists == TRUE
              ) %>% 
              group_by(
                MF,
                Adduct
              ) %>% 
              arrange(status) %>% 
              slice(1)
            
            statuses %>% 
              rowwise() %>% 
              group_split() %>%
              map_dfr(
                ~fileName(.x$MF,.x$database,.x$organism,path = library_path) %>% 
                  read_rds() %>% 
                  consensus(
                    adduct = .x$Adduct[1],
                    threshold = threshold
                  ) %>% 
                  mutate(
                    mf = .x$MF
                  ) %>% 
                  relocate(
                    mf,
                    .before = adduct
                  )
              )
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
