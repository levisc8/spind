wrm.pad<-function(padform=0,padzone=1){
  ###############################################################################
  # Description
  # Auxiliary function for wavelet coefficient padding. 
  # Typically only used internally by WRM.
  # Arguments: 
  # padform   0 for padding with zeros, 1 for padding with mean values, 
  #           2 for padding with mirror values. Padform is automatically set to 
  #           zero in case of either level=0 or a formula including an intercept 
  #           and a non-gaussian family
  # padzone   factor for expanding the padding zone
  #
  # Value:     A list with components named as the arguments. 
  ###############################################################################
  pad.list<-list(padform=padform,padzone=padzone)
  return(pad.list)
}