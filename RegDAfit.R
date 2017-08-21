#Fitting a RegDA model
#Input:
  #Input: covarites
  #Response: response variable
  #alpha: 1st order f regularization (balance between LDA & QDA models)
  #gamma: 2nd order of regularization (shrinking generated covariance matrix toward scalar diagonal)
#Output:
  #yhat:prediction
  #Mu_k: means matrix
  #sigmaHatMin1_k: inverse of each class covariance matrix
  #Pi_k=a priori probs
  #det_k=determinant of each class covariate matrix

RegDAfit<- function(Input,Response,alpha,gamma){
  Input=as.matrix(Input)
  N=nrow(Input)
  p=ncol(Input)
  K=length(unique(Response))
  classes=as.numeric(levels(as.factor(Response)))
  #init parameters
  Mu_k=matrix(0,nrow = p,ncol = K)#centroid vectors
  Pi_k=vector(length= K)# classes a priori
  sigma=matrix(0,nrow=p,ncol=p)# within class covariance (constant)
  sigma_k=list()#within class covariance
  sigmaHat=matrix(0,ncol=p,nrow = p)
  det_k=vector(length = K)
  #compute parameters
  for (i in 1:K){
    sigma_k[[i]]=matrix(0,ncol=p,nrow = p)
    tmp=Input[Response==classes[i],]
    N0=nrow(tmp)
    Pi_k[i]=N0/N
    Mu_k[,i]=colMeans(tmp)
    for(j in 1:N0){
      sigma_k[[i]]=sigma_k[[i]]+(tmp[j,]-Mu_k[,i])%*%t(tmp[j,]-Mu_k[,i])
      sigma=sigma+(tmp[j,]-Mu_k[,i])%*%t(tmp[j,]-Mu_k[,i])
    }
    sigma_k[[i]]=(1/(N0-1))*sigma_k[[i]]
  }
  sigma=(1/(N-K))*sigma
  sigma2=(1/p)*sum(diag(sigma))
  sigmaHat=gamma*sigma+(1-gamma)*sigma2*diag(p)
  sigmaHat_k=list()
  sigmaHatMin1_k=list()
  for (i in 1:K){
    sigmaHat_k[[i]]=matrix(0,ncol=p,nrow=p)
    sigmaHatMin1_k[[i]]=matrix(0,ncol=p,nrow=p)
    sigmaHat_k[[i]]=alpha*sigma_k[[i]]+(1-alpha)*sigmaHat
    sigmaHatMin1_k[[i]]=solve(sigmaHat_k[[i]])
    det_k[i]=det(sigmaHat_k[[i]])
  }
  #Classification
  deltaTrain=matrix(0,ncol = K,nrow = N)
  for (i in 1:K){
    deltaTrain[,i]=apply(Input,1,function(x) -0.5*t(x-Mu_k[,i])%*%sigmaHatMin1_k[[i]]%*%(x-Mu_k[,i]))
    deltaTrain[,i]=deltaTrain[,i]-0.5*log(det_k[i])+log(Pi_k[i])
  }
  #prediction
  yHatTrain=apply(deltaTrain,1,which.max)
  return(list(yhat=yHatTrain,Mu_k=Mu_k,sigmaHatMin1_k=sigmaHatMin1_k,Pi_k=Pi_k,det_k=det_k))
}
