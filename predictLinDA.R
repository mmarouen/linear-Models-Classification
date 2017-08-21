#prediction function for LDA
#Input:
  #Input: covariates
  #modlinda:output of LDA model
#output:
  #predictions

predictLinDA<-function(Input,modLinda){
  sigmaMinus1=modLinda$sigmaMin1
  centroidMat=modLinda$centroids
  Pi_k=modLinda$Pi_k
  K=ncol(centroidMat)
  #compute discriminant matrix
  deltaMat=matrix(0,ncol=K,nrow=nrow(Input))
  for(i in 1:K){
    deltaMat[,i]=as.matrix(Input)%*%sigmaMinus1%*%centroidMat[,i]
    deltaMat[,i]=deltaMat[,i]-0.5*t(centroidMat[,i])%*%sigmaMinus1%*%centroidMat[,i]+log(Pi_k[i])
  }
  yHat=apply(deltaMat,1, which.max)
  return(yHat)
}
