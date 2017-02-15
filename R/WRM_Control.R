wrm.control<-function(eps=1e-5,denom.eps=1e-20,itmax=200){
  ###############################################################################
  # Description
  # Auxiliary function for WRM fitting. Typically only used internally by WRM.
  # Arguments:
  # eps       positive convergence tolerance
  # denom.eps the iterations converge when
  #           max(abs(coeff.new-coeff.old)/(abs(coeff.old)+denom.eps)) <= eps
  # itmax     integer giving the maximal number of iterations
  #
  # Value:     A list with components named as the arguments.
  ###############################################################################
  contr.list<-list(eps=eps,denom.eps=denom.eps,itmax=itmax)
  return(contr.list)
}
#
