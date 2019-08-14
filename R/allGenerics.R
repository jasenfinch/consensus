#' @rdname keggConsensus

setGeneric('keggConsensus',function(x,organism = 'hsa', threshold = 0.5){
  standardGeneric('keggConsensus')
})

#' @rdname construct

setGeneric('construct',function(x,organism = 'hsa', threshold = 0.5, databases = c('kegg','pubchem')){
  standardGeneric('construct')
})