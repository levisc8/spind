wrm.moran<-function(lim1=0,increment=1){
  ###############################################################################
  # Description
  # Auxiliary function for specifying the Moran's I calculation. 
  # Typically only used internally by WRM.
  # Arguments: 
  # lim1      lower limit for first bin
  # increment increment (=1 and lim1=0 is conform to correlog{ncf})
  #
  # Value:     A list with components named as the arguments. 
  ###############################################################################
  moran.list<-list(lim1=lim1,increment=increment)
  return(moran.list)
}