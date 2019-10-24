#' mf
#' @rdname mf
#' @description Get and the molecular formula of a Consensus object. 
#' @param x S4 object of class Consensus
#' @export

setMethod('mf',signature = 'Consensus',
          function(x){
            x@MF
          })

#' adductRules
#' @rdname adductRules
#' @description Get the adduct rules of a Consensus object. 
#' @param x S4 object of class Consensus
#' @export

setMethod('adductRules',signature = 'Consensus',
          function(x){
            x@adductRules
          })

#' organism
#' @rdname organism
#' @description Get the organism ID of a Consensus object. 
#' @param x S4 object of class Consensus
#' @export

setMethod('organism',signature = 'Consensus',
          function(x){
            x@organism
          })

#' database
#' @rdname database
#' @description Get the database of a Consensus object. 
#' @param x S4 object of class Consensus
#' @export

setMethod('database',signature = 'Consensus',
          function(x){
            x@database
          })

#' threshold
#' @rdname threshold
#' @description Get the threshold of a Consensus object. 
#' @param x S4 object of class Consensus
#' @export

setMethod('threshold',signature = 'Consensus',
          function(x){
            x@threshold
          })

#' hits
#' @rdname hits
#' @description Get the database molecular formula matches of a Consensus object. 
#' @param x S4 object of class Consensus
#' @export

setMethod('hits',signature = 'Consensus',
          function(x){
            x@hits
          })

#' PIPs
#' @rdname PIPs
#' @description Get the putative ionisation products of a Consensus object. 
#' @param x S4 object of class Consensus
#' @export

setMethod('PIPs',signature = 'Consensus',
          function(x){
            x@PIPs
          })

#' classifications
#' @rdname classifications
#' @description Get the classifications of a Consensus object. 
#' @param x S4 object of class Consensus
#' @export

setMethod('classifications',signature = 'Consensus',
          function(x){
            x@classifications
          })

#' consensusClassifications
#' @rdname consensusClassifications
#' @description Get the consensus classifications of a Consensus object. 
#' @param x S4 object of class Consensus
#' @export

setMethod('consensusClassifications',signature = 'Consensus',
          function(x){
            x@consensus
          })