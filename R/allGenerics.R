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

setGeneric('mfHits',function(x){
  standardGeneric('mfHits')
})

setGeneric('pips',function(x){
  standardGeneric('pips')
})

setGeneric('classify',function(x){
  standardGeneric('classify')
})

setGeneric('consensus',function(x){
  standardGeneric('consensus')
})

setGeneric('overallConsensus',function(x){
  standardGeneric('overallConsensus')
})

setGeneric('classificationTree',function(x){
  standardGeneric('classificationTree')
})

setGeneric('saveConsensus',function(x,path = '.'){
  standardGeneric('saveConsensus')
})