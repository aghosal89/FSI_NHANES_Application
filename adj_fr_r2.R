## This function computes the Frechet adjusted R^2 for predictions for quantiles 

#  INPUTS:
#
#   qin     -  an nxm matrix whose rows are the quantiles corresponding to a distribution
#   qpred   -  an nxm matrix whose rows are predicted quantiles for the Y above. 
#   tt      -  an equidistant grid of points on [0,1] of length m for quantiles.
#   q       -  number of covariates in the model. 
#
#  OUTPUT:  -  A list contanining Frechet R-squared, Adjusted Frechet R-squared.

adj_fr_r2 <- function(qin=NULL, qpred=NULL, tt=NULL, q=NULL) {
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
  frmean<- apply(qpred, 2, mean)
  
  library("fdadensity")
  ss1<- mean(sapply(1:n, function(i) fdadensity:::trapzRcpp(X=tt, Y=(qin[i,] - qpred[i,])^2)))
  ss2<- mean(sapply(1:n, function(i) fdadensity:::trapzRcpp(X=tt, Y=(qin[i,] - frmean)^2)))
  fr_r2 <- 1-(ss1/ss2)
  adj_fr_r2 <- fr_r2 - (1 - fr_r2)*(q/(n-q-1))
  return(list(Frechet_R2=fr_r2,Frechet_Adj_R2=adj_fr_r2))
}
