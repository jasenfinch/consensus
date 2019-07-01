
setGeneric('pubchemPIPs',function(ips){
  standardGeneric('pubchemPIPs')
})

setGeneric('keggConsensus',function(x,organism = 'hsa'){
  standardGeneric('keggConsensus')
})

setGeneric('consensus',function(x,filterUnclassified = F){
  standardGeneric('consensus')
})