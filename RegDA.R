#Regularized discriminant analysis model
#Input:
  #Input:covariates
  #Response: response vector
  #alpha: 1st regularization parameter (balance between an LDA model and QDA model)
  #gamma: 2nd regularization parameter (shrinking generated covariance toward scalar diagonal)
  #CV: binary for CV
  #cvScale:binary for scaling during cv
  #L:resolution of regularization

RegDA<-function(Input,Response,alpha=NULL,gamma=NULL,CV=FALSE,cvScale=FALSE,L=40){
  
  #init parameters:
  error=c()
  errorMat=c()#selection error matrix
  yhat=c()#training prediction
  stderr=c()#standard deviations for CV

  if(!is.null(alpha) & !is.null(gamma)){
    out1=RegDAfit(Input,Response,alpha=alpha,gamma=gamma)
    i0=1
    j0=1
  }
  if(is.null(alpha) | is.null(gamma)){
    if(is.null(alpha)){alpha=seq(0,1,length=L)}
    if(is.null(gamma)){gamma=seq(0,1,length=L)}
    errorMat=matrix(0,ncol=length(gamma),nrow=length(alpha))
    stderr=matrix(0,ncol=length(gamma),nrow=length(alpha))
    for(i in 1:length(alpha)){
      for(j in 1:length(gamma)){
        if(!CV){
          out=RegDAfit(Input,Response,alpha=alpha[i],gamma=gamma[j])
          errorMat[i,j]=mean(Response != predictRegDA(Input,out))
        }
        if(CV){
          source("E:/R Project/toolkits/ModelAssessment&Selection.R")
          out=crossValidate(Input,Response,type = "RegDA",complexity = c(alpha[i],gamma[j]),scaleInputs = cvScale,scaleType = "Standardize")
          errorMat[i,j]=mean(out$errorVector)
          stderr[i,j]=out$sdVal
        }
      }
    }
    minInd=which(error==min(error),arr.ind = TRUE)
    i0=minInd[1]
    j0=minInd[2]
    out1=RegDAfit(Input,Response,alpha = alpha[i0],gamma=gamma[j0])
  }
  yhat=predictRegDA(Input,out1)
  error=mean(yhat != Response)
  return(list(yhat=yhat,error=error,alphaVal=alpha[i0],gammaVal=gamma[j0],modRegda=out1))
}
