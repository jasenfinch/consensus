#' @rdname keggConstruct

setGeneric('keggConstruct',function(x,organism = 'hsa', threshold = 50){
  standardGeneric('keggConstruct')
})

#' @rdname construct

setGeneric('construct',function(x,organism = 'hsa', threshold = 50, databases = c('kegg','pubchem')){
  standardGeneric('construct')
})