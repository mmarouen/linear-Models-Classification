#fitting a given RRDA model
#Input:
  #Input:covariates
  #Response: response vector
  #sub: order of regularization (size of covariance matrix)
  #modelF:if TRUE then a best model is returned
#Output:
  #yHatTrain:prediction
  #error:training error
  #sub:best order
  #canonicalMeans: canonical means matrix (in projected space)
  #canonicalMat: canonical covariance matrix (in projected space)
  #Pi_k: means
  #canonicalTrainData: canonical data matrix

RRDAfit <- function(Input,Response,sub=NULL,modelF=FALSE){
  
  Input=as.matrix(Input)
  N=nrow(Input)
  p=ncol(Input)
  K=length(unique(Response))
  classes=sort(unique(yTrain))
  error=c()
  yHatTrainMat=matrix()
  yHatTrain=c()

  if(is.null(sub)){modelF=TRUE}
  #initiate variables
  Mu_k=matrix(0,nrow = K,ncol = p)#centroid vectors
  Pi_k=vector(length= K)# classes a priori
  sigma=matrix(0,nrow=p,ncol=p)# within class common covariance 
  mu=colMeans(Mu_k)#general mean
  B=matrix(0,ncol=p,nrow = p) #between class covariance B
  
  #compute variables in features space
  for (i in 1:K){
    tmp=Input[Response==classes[i],]
    N0=nrow(tmp)
    Pi_k[i]=N0/N
    Mu_k[i,]=colMeans(tmp)
    for(j in 1:N0){
      sigma=sigma+(tmp[j,]-Mu_k[i,])%*%t(tmp[j,]-Mu_k[i,])}
  }
  sigma=(1/(N-K))*sigma#within class covariance
  for(i in 1:nrow(Mu_k)){
    B=B+Pi_k[i]*(Mu_k[i,]-mu)%*%t(Mu_k[i,]-mu)
  }#between class covariance
  
  #Sphere data w.r.t sigma 
  e=eigen(sigma)#sigma=V*D*t(V)
  U=e$vectors
  D=diag(e$values)
  W_minus_half=U%*%diag(sqrt(e$values)^(-1))
  Mstar=Mu_k%*%W_minus_half
  Inputstar=Input%*%W_minus_half
  Bstar=t(W_minus_half)%*%B%*%W_minus_half
  
  #extract eigen vectors of Bstar
  e2=eigen(Bstar)
  V=e2$vector
  
  #Project into canonical base
  InputTrainCan=Inputstar%*%V
  MeanCan=Mstar%*%V
  canonicalMatrix=W_minus_half%*%V
  #train classifier
  deltaTrain=list()

  if(modelF){
    yHatTrainMat=matrix(nrow = nrow(Input),ncol = p)
    for(subspace in 1:p){
      deltaTrain[[subspace]]=matrix(0,ncol=K,nrow=N)
      for (i in 1:K){
        for (m in 1:N){
          deltaTrain[[subspace]][m,i]=0.5*sum((InputTrainCan[m,1:subspace]-MeanCan[i,1:subspace])^2)-log(Pi_k[i])
        }
      }
      yHatTrainMat[,subspace]=apply(deltaTrain[[subspace]],1,which.min)
      error[subspace]=mean(yHatTrainMat[,subspace] != Response)
    }
    sub=which.min(error)
    yHatTrain=yHatTrainMat[,sub]
  }
  
  if(!modelF){
    deltaTrain[[1]]=matrix(0,ncol=K,nrow=N)
    for (i in 1:K){
      for (m in 1:N){
        deltaTrain[[1]][m,i]=0.5*sum((InputTrainCan[m,1:sub]-MeanCan[i,1:sub])^2)-log(Pi_k[i])
      }
    }
    yHatTrain=apply(deltaTrain[[1]],1,which.min)
  }
  return(list(yHatTrain=yHatTrain,error=error,sub=sub,canonicalMeans=MeanCan,canonicalMat=canonicalMatrix,Pi_k=Pi_k,canonicalTrainData=InputTrainCan))
}
