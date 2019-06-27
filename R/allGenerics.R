
setGeneric('pubchemPIPs',function(ips){
  standardGeneric('pubchemPIPs')
})

setGeneric('pipClassifications',function(ips,nCores = detectCores(), clusterType = 'PSOCK'){
  standardGeneric('pipClassifications')
})