#fits a RRDA model
#Input:
  #Input: covariates
  #Response: response vector
  #sub:regularization order
  #CV: binary for CV
  #cvScale: binary for regularizing during CV
#Output:
  #yhat: prediction
  #error: training error
  #errorMat: error matrix in case of CV
  #stderr: errors stddev in case CV
  #sub: optimal regularization order
  #canonicalMatrix: canonical covariance matrix
  #Pi_k: a priori probs
  #canonicalMeans: canonical means matrix

RRDA<-function(Input,Response,sub=NULL,CV=FALSE,cvScale=FALSE){
  
  #init parameters:
  error=c()
  errorMat=c()#selection error matrix
  yhat=c()#training prediction
  stderr=c()#standard deviations for CV
  if(!is.null(sub)){
    out=RRDAfit(Input,Response,sub = sub,modelF=FALSE)
  }
  if(is.null(sub)){
    if(!CV){
      out0=RRDAfit(Input,Response,modelF=TRUE)
      sub=out0$sub
      out=RRDAfit(Input,Response,sub=sub,modelF=FALSE)
    }
    if(CV){
      source("D:/RProject/toolkits/ModelAssessment&Selection.R")
      p=ncol(Input)
      for(i in 1:p){
        out0=crossValidate(Input,Response,type="RRDA",complexity = i,scaleInputs = cvScale,scaleType = "Standardize")
        errorMat[i]=mean(out$errorVector)
        stderr[i]=out$sdVal
      }
      sub=which.min(error)
      out=RRDAfit(Input,Response,sub=sub,modelF=FALSE)
    }
  }
  canonicalMeans=out$canonicalMeans
  Pi_k=out$Pi_k
  canonicalMatrix=out$canonicalMat
  sub=out$sub
  yhat=predictRRDA(Input,out)
  error=mean(yhat != Response)
  return(list(yhat=yhat,error=error,errorMat=errorMat,stderr=stderr,sub=sub,canonicalMatrix=canonicalMatrix,Pi_k=Pi_k,canonicalMeans=canonicalMeans))
}
