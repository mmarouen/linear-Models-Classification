#linear classifier
#Input:
  #Input:covariates
  #model: output of a linear classification model
#Output: classification

lcPredict<-function(Input,model){
  fit=model$fit
  if(fit=="LDA") return(predictLinDA(Input,model))
  if(fit=="QDA") return(predictQuaDA(Input,model))
  if(fit=="RRDA")return(predictRRDA(Input,model))
  if(fit=="RegDA") return(predictRegDA(Input,model))
}
