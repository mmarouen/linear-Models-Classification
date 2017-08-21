#prediction for RRDA model
#Input:
  #Input: covariates
  #modRRDA: RRDA model
#Output:yhat=prediction

predictRRDA<-function(Input,modRRDA){

  #load model parameters
  sub=modRRDA$sub
  meancan=modRRDA$canonicalMeans
  Pi_k=modRRDA$Pi_k
  canonicalMatrix=as.matrix(modRRDA$canonicalMat)
  K=length(Pi_k)
  #discriminant matrix
  N=nrow(Input)
  canonicalInput=as.matrix(Input)%*%canonicalMatrix
  discrMat=matrix(0,ncol=K,nrow = N)
  
  for (i in 1:K){
    for (m in 1:N){
      discrMat[m,i]=0.5*sum((canonicalInput[m,1:sub]-meancan[i,1:sub])^2)-log(Pi_k[i])
    }
  }
  yHat=apply(discrMat,1,which.min)
  return(yHat)
}
