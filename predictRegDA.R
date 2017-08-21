#prediction for RegDA model
#Input:
  #Input:covariates
  #modelRegDA: output of RegDA model
#Output: yhat=prediction

predictRegDA<-function(Input,modeRegDA){
  #load model parameters
  if(!is.null(modeRegDA$modRegda)){modeRegDA=modeRegDA$modRegda}
  sigmaMinus1_k=modeRegDA$sigmaHatMin1_k
  Mu_k=modeRegDA$Mu_k
  det_k=modeRegDA$det_k
  Pi_k=modeRegDA$Pi_k
  K=length(Pi_k)
  #discriminant matrix
  descrMat=matrix(0,ncol=K,nrow=nrow(Input))
  for(i in 1:K){
    descrMat[,i]=apply(as.matrix(Input),1,function(x) -0.5*t(x-Mu_k[,i])%*%sigmaMinus1_k[[i]]%*%(x-Mu_k[,i]))
    descrMat[,i]=descrMat[,i]-0.5*log(det_k[i])+log(Pi_k[i])
  }
  yHat=apply(descrMat,1,which.max)
  return(yHat)
}
