#Quadratic discriminant analysis
#Input:
  #Input: covariates
  #Response: response vector
#Output:
  #yhat:prediction
  #sigmaMinus1_k: covariance matrices inverses
  #Mu_k: centroids matrix
  #det_k: coavriance matrices determinant

QuaDA<- function(Input,Response){
  Input=as.matrix(Input)
  N=nrow(Input)
  p=ncol(Input)
  K=length(unique(Response))
  classes=as.numeric(levels(as.factor(Response)))
  Mu_k=matrix(0,nrow = p,ncol = K)#centroid vectors
  Pi_k=vector(length= K)# classes a priori
  sigma_k=list()#data covariances matrix per class
  sigmaMinus1_k=list()
  det_k=vector(length = K)
  for (i in 1:K){
    sigma_k[[i]]=matrix(0,ncol=p,nrow = p)
    tmp=Input[Response==classes[i],]
    N0=nrow(tmp)
    Pi_k[i]=N0/N
    Mu_k[,i]=colMeans(tmp)
    for(j in 1:N0){
      sigma_k[[i]]=sigma_k[[i]]+(tmp[j,]-Mu_k[,i])%*%t(tmp[j,]-Mu_k[,i])
    }
    sigma_k[[i]]=(1/(N0-1))*sigma_k[[i]]
    sigmaMinus1_k[[i]]=solve(sigma_k[[i]])
    det_k[i]=det(sigma_k[[i]])
  }
  #discriminant matrix
  deltaTrain=matrix(0,ncol = K,nrow = N)
  for (i in 1:K){
    deltaTrain[,i]=apply(Input,1,function(x) -0.5*t(x-Mu_k[,i])%*%sigmaMinus1_k[[i]]%*%(x-Mu_k[,i]))
    deltaTrain[,i]=deltaTrain[,i]-0.5*log(det_k[i])+log(Pi_k[i])
  }
  #prediction
  yHatTrain=apply(deltaTrain,1,which.max)
  return(list(yhat=yHatTrain,sigmaMinus1_k=sigmaMinus1_k,Mu_k=Mu_k,det_k=det_k,Pi_k=Pi_k))
}
