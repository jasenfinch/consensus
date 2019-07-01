
setGeneric('pubchemPIPs',function(ips){
  standardGeneric('pubchemPIPs')
})

setGeneric('keggConsensus',function(x,organism = 'hsa', threshold = 0.5){
  standardGeneric('keggConsensus')
})

setGeneric('consensus',function(x,organism = 'hsa', threshold = 0.5){
  standardGeneric('consensus')
})