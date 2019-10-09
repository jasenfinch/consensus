#' @rdname mf
setGeneric('mf',function(x){
  standardGeneric('mf')
})

#' @rdname adductRules
setGeneric('adductRules',function(x){
  standardGeneric('adductRules')
})

#' @rdname organism
setGeneric('organism',function(x){
  standardGeneric('organism')
})

#' @rdname database
setGeneric('database',function(x){
  standardGeneric('database')
})

#' @rdname threshold
setGeneric('threshold',function(x){
  standardGeneric('threshold')
})

#' @rdname hits
setGeneric('hits',function(x){
  standardGeneric('hits')
})

#' @rdname PIPs
setGeneric('PIPs',function(x){
  standardGeneric('PIPs')
})

#' @rdname classifications
setGeneric('classifications',function(x){
  standardGeneric('classifications')
})

#' @rdname consensusClassifications
setGeneric('consensusClassifications',function(x){
  standardGeneric('consensusClassifications')
})

#' @rdname mfHits
setGeneric('mfHits',function(x){
  standardGeneric('mfHits')
})

#' @rdname keggConstruct

setGeneric('keggConstruct',function(x,organism = 'hsa', threshold = 50){
  standardGeneric('keggConstruct')
})

