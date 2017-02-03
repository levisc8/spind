padding<-function(M){
  ###############################################################################
  # Description
  # A function for mirror padding of a matrix 
  # Arguments:
  # M       a square matrix
  # Value:  padding returns a new matrix
  ###############################################################################
  if(dim(M)[1]!=dim(M)[2]) stop("no square matrix")
  n<-dim(M)[1]
  P<-which(!is.na(M), arr.ind = TRUE)
  # columns
  for (i in 1:n){
    if(length(intersect(P[,2],i))==0) next
    mi<-min(P[which(P[,2]==i),1]) 
    ma<-max(P[which(P[,2]==i),1])
    lmi<-mi-1
    lmis<-mi:(mi+lmi-1)
    if(mi>1 & length(lmis[lmis<n])==length(lmis)) M[lmi:1,i]<-M[mi:(mi+lmi-1),i]
    if(mi>1 & length(lmis[lmis<n])!=length(lmis)) 
      M[lmi:(lmi-length(mi:n)+1),i]<-M[mi:n,i]
    lma<-ma+1
    lmas<- ma:(ma+lma-n)
    if(ma<n & length(lmas[lmas>0])==length(lmas)) M[lma:n,i]<-M[ma:(ma+lma-n),i]
    if(ma<n & length(lmas[lmas>0])!=length(lmas)) 
      M[lma:(lma+length(ma:1)-1),i]<-M[ma:1,i]
  }
  # transpose for rows
  M<-t(M)
  P<-which(!is.na(M), arr.ind = TRUE)
  # rows 
  for (i in 1:n){
    if(length(intersect(P[,2],i))==0) next
    mi<-min(P[which(P[,2]==i),1]) 
    ma<-max(P[which(P[,2]==i),1])
    lmi<-mi-1
    lmis<-mi:(mi+lmi-1)
    if(mi>1 & length(lmis[lmis<n])==length(lmis)) M[lmi:1,i]<-M[mi:(mi+lmi-1),i]
    if(mi>1 & length(lmis[lmis<n])!=length(lmis)) 
      M[lmi:(lmi-length(mi:n)+1),i]<-M[mi:n,i]
    lma<-ma+1
    lmas<- ma:(ma+lma-n)
    if(ma<n & length(lmas[lmas>0])==length(lmas)) M[lma:n,i]<-M[ma:(ma+lma-n),i]
    if(ma<n & length(lmas[lmas>0])!=length(lmas)) 
      M[lma:(lma+length(ma:1)-1),i]<-M[ma:1,i]
  }
  # back transform
  M<-t(M)
  M[is.na(M)]<-mean(M, na.rm=TRUE)
  M
}