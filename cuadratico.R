## This function projects the model predictions into the L^2-Wasserstein space of distributions

## Inputs:  

# 1) prediciones  - model predictions of the GLM
# 2) cotainferior - lower limit of the quantiles
# 3) cotasuperior - upper limit of the quantiles

## Output:  

# a matrix of same dimension as the input matrix prediciones, of which the rows are projected onto the space of quantiles.

## R libraries to be read to run the function

# 1) osqp
# 2) Rcpp
# 3) Matrix

cuadratico= function(prediciones, cotainferior=-10e-5, cotasuperior=800) {
  prediciones= as.matrix(prediciones)
  n= dim(prediciones)[1]
  p= dim(prediciones)[2]
  
  salida= matrix(0, nrow= n, ncol= p)
  
  for(i in 1:n){
    
    P= diag(p)
    A= diag(p)*-1
    for(j in 1:(p-1)) {
      A[j,j+1]=1
    }
    A[p,p]=0
    u= rep(cotainferior,p-1)
    l= rep(cotasuperior,p-1)
    u= c(u,cotainferior)
    l= c(l,cotasuperior)
    q= -prediciones[i,]
    # u, l change
    res <- do.call(osqp::solve_osqp, list(P = P, q = q, 
                  A = A, l =u, u=l, pars = osqp::osqpSettings(verbose = FALSE)))
    res= sort(res$x)
    salida[i,]= res
  }
  return(salida)
}

