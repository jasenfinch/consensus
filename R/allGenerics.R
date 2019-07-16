#' @rdname keggConsensus

setGeneric('keggConsensus',function(x,organism = 'hsa', threshold = 0.5){
  standardGeneric('keggConsensus')
})

#' @rdname consensus

setGeneric('consensus',function(x,organism = 'hsa', threshold = 0.5){
  standardGeneric('consensus')
})