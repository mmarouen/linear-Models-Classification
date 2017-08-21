#linear classification function supports: LDA, QDA, RRDA, RDA
#Input:
  #Input: input covariates
  #Response: response variable
  #fit: model to be applied
  #sub,alpha,gamma: different regularization parameters for models
  #cv: binary for cross validation
  #cvScale: binary, for scaling variables prior to cv or not
#Output:
  #results specific to model
  #fit:string containing model name
lc<-function(Input,Response,fit="LDA",sub=NULL,alpha=NULL,gamma=NULL,cv=FALSE,cvScale=FALSE){

  ll=list()
  if(fit=="LDA"){
    ll=LinDA(Input,Response)
  }
  if(fit=="QDA"){
    ll=QuaDA(Input,Response)
  }
  if(fit=="RRDA"){
    ll=RRDA(Input,Response,sub=sub,CV = cv,cvScale = cvScale)
  }
  if(fit=="RegDA"){
    ll=RegDA(Input,Response,alpha = alpha,gamma = gamma,CV = cv,cvScale = cvScale)
  }
  ll$fit=fit
  return(ll)
}
