#Prediction for QDA model
#Input:
  #Input: covariates matrix
  #modQuaDA: QDA model
#Output: yhat:prediction

predictQuaDA<-function(Input,modQuaDA){
  #load model parameters
  sigmaMinus1_k=modQuaDA$sigmaMinus1_k
  Mu_k=modQuaDA$Mu_k
  det_k=modQuaDA$det_k
  Pi_k=modQuaDA$Pi_k
  K=length(Pi_k)
  #discriminant matrix
  descrMat=matrix(0,ncol=K,nrow=nrow(Input))
  for(i in 1:K){
    descrMat[,i]=apply(Input,1,function(x) -0.5*t(x-Mu_k[,i])%*%sigmaMinus1_k[[i]]%*%(x-Mu_k[,i]))
    descrMat[,i]=descrMat[,i]-0.5*log(det_k[i])+log(Pi_k[i])
  }
  yHat=apply(descrMat,1,which.max)
  return(yHat)
}
