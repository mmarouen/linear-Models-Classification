#Linear Discriminant Analysis model
#Input:
  #Input: covariates
  #Response: response vector
#Output:
   #yhat: prediction
   #sigmaMin1: inverse of covariance matrix
   #centroids: centroids matrix
   #Pi_k: a priori probs

LinDA<- function(Input,Response){
  Input=as.matrix(Input)
  N=nrow(Input)
  K=length(unique(Response))
  classes=as.numeric(levels(as.factor(Response)))
  p=ncol(Input)
  Mu_k=matrix(0,nrow = p,ncol = K)#centroid vectors
  Pi_k=vector(length= K)# a priori probabilities
  sigma=matrix(0,ncol=p,nrow = p)#classes covariance matrix
  for (i in 1:K){
    tmp=Input[Response==classes[i],]
    Pi_k[i]=nrow(tmp)/N
    Mu_k[,i]=colMeans(tmp)
    for(j in 1:nrow(tmp)){
      sigma=sigma+(tmp[j,]-Mu_k[,i])%*%t(tmp[j,]-Mu_k[,i])
    }
  }
  sigma=sigma*(1/(N-K))
  sigmaMinus1=solve(sigma)
  
  #Classification
  deltaTrain=matrix(0,ncol = K,nrow = N)
  for (i in 1:K){
    deltaTrain[,i]=Input%*%sigmaMinus1%*%Mu_k[,i]
    deltaTrain[,i]=deltaTrain[,i]-0.5*t(Mu_k[,i])%*%sigmaMinus1%*%Mu_k[,i]+log(Pi_k[i])
  }
  #prediction
  yHatTrain=apply(deltaTrain,1, which.max)
  return(list(yhat=yHatTrain,sigmaMin1=sigmaMinus1,centroids=Mu_k,Pi_k=Pi_k))
}
