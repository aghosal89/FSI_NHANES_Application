## This function computes the Fréchet R^2 and the adjusted Fréchet R^2 for predicticted quantiles of 
# the physical activity distributions.

# R library required to run this function:
# fdadensity

#  INPUTS:
#
#   qin             -  an nxm matrix whose rows are the quantiles corresponding to a distribution
#   qpred           -  an nxm matrix whose rows are predicted quantiles for the Y above. 
#   tt              -  an equidistant grid of points on [0,1] of length m for quantiles.
#   q               -  number of covariates in the model. 
#   survey_weights  -  Survey weights used for the individuals in NHANES study
#
#  OUTPUT  -  A list contanining Fréchet R-squared (Frechet_R2), 
#                     Adjusted Fréchet R-squared (Adj_Frechet_R2).

adj_fr_r2 <- function(qin=NULL, qpred=NULL, tt=NULL, q=NULL, survey_weights= NULL) {
  
  # perform checks
  if(nrow(qin)!=nrow(qpred) | ncol(qin)!=ncol(qpred)) {
    stop('Dimentions of response and its predictions do not match.')
  }
  if(!is.matrix(qin) | !is.matrix(qpred)) {
    message('the response and its predictions should be matrices')
    qin <- as.matrix(qin)
    qpred <- as.matrix(qpred)
  }
  n <- nrow(qin)
  
  w<- survey_weights/sum(survey_weights)
  frmean<- apply(w*qpred, 2, sum) 
  
  ss1<- sapply(1:n, function(i) fdadensity:::trapzRcpp(X=tt, Y=(qin[i,] - qpred[i,])^2))
  ss2<- sapply(1:n, function(i) fdadensity:::trapzRcpp(X=tt, Y=(qin[i,] - frmean)^2))
  
  fr_r2 <- 1-(sum(w*ss1)/sum(w*ss2))
  adj_fr_r2 <- fr_r2 - (1 - fr_r2)*(q/(n-q-1))
               
  return(list(Frechet_R2=fr_r2, Adj_Frechet_R2=adj_fr_r2))
  
}
