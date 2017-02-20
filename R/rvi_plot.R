# Model selection tables and relative variable importance
# in terms of scale levels
#
# Description The summary is based on the results provided by
# multi-model inference for wavelet multiresolution regression (mmiWMRR).
# Based on the results of this approach, the relative
# importance of predictors are calculated. The scale-dependent results
# are graphically displayed.
#
# Arguments
# formula   ...
# family
# data
# coord
# maxlevel   an integer for maximum scale level
# detail
# wavelet
# wtrafo
# n.eff


rvi.plot<-function(formula,family,data,coord,maxlevel,detail=TRUE,
wavelet="haar",wtrafo="dwt",n.eff=NULL){

cat("\n","Model selection tables:","\n","\n")

mmi<- mmiWMRR(formula,family,data,coord,scale=1,
detail=detail,wavelet=wavelet,wtrafo=wtrafo,n.eff=n.eff)
nrowA<-dim(mmi$result)[1]
ncolA<-dim(mmi$result)[2]

nvar<-dim(mmi$result)[2]-6
leg<-dimnames(mmi$result)[[2]][2:(nvar+1)]

A<-array(NA,c(nrowA,ncolA,maxlevel))
level<-rep(NA,maxlevel)
A[,,1]<-mmi$result
level[1]<-mmi$level

if(maxlevel>=2){
for (i in 2:maxlevel) {
mmi<- mmiWMRR(formula,family,data,coord,scale=i,
detail=detail,wavelet=wavelet,wtrafo=wtrafo,n.eff=n.eff)
A[,,i]<-mmi$result
level[i]<-mmi$level
}}

# Plot: scale-dependent relative variable importance
cat("\n","---","\n","Relative variable importance:","\n","\n")

klimitscale<-dim(A)[3]
ip<-dim(A)[1]

WeightSums<-matrix(NA,nvar,klimitscale)
for (kscale in 1:klimitscale){
for(kvar in 2:(nvar+1)){
for (i in 1: ip){
if(!is.na(A[i,kvar,kscale])) A[i,kvar,kscale]<-A[i,(nvar+6),kscale]
}}
B<-A[1:ip,2:(nvar+1),kscale]
WeightSums[,kscale]<-colSums(B,na.rm=TRUE)
} # kscale
vec<-1:nvar
plot(level,WeightSums[1,],type="b",ylim=c(0,2),
xlim=c(min(level),max(level)),ylab="Relative Variable Importance",
xlab="Level", pch=2,lty=vec[1],lwd=2)
for (kvar in 2:nvar) {
points(level,WeightSums[kvar,],type="b",pch=kvar+1,
lty=vec[kvar],lwd=2)
}
#leg<-dimnames(mmi$res)[[2]][vec+1]
v<-2:(nvar+1)
legend(2,2,leg,pch=v,lty=vec,lwd=2)

rownames(WeightSums)<-leg
colnames(WeightSums)<-paste("level",c(1:klimitscale),sep="=")
print(WeightSums)

fit<-list(rvi=WeightSums)
fit

}


