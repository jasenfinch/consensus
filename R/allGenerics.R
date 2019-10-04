#' @rdname keggConstruct

setGeneric('keggConstruct',function(x,organism = 'hsa', threshold = 0.5){
  standardGeneric('keggConstruct')
})

#' @rdname construct

setGeneric('construct',function(x,organism = 'hsa', threshold = 0.5, databases = c('kegg','pubchem')){
  standardGeneric('construct')
})