
# This function fits the partially-linear Frechet Single Index model to the 
# distributional responses. The files 'wn_cost.R', 'polar2cart.R', 'cart2polar.R',
# 'survey2wassersteinmodel.R', 'cuadratico.R' have to be sourced prior to running 
# the codes in this script.

## Inputs:
# 
# tt          - length m grid spanning [0, 1], used as grid for quantile functions
# datosfda    - nxm matrix of response quantile functions on grid tt
# si_vars     - A p-vector of variables' names to be considered in the Single Index part.
# linear_vars - A q-vector of variables' names to be considered in the linear part.
# formula_lv  - a single character denoting the formula for the variables in the 
#               linear part to be considered for the model, e.g. interactions, 
#               higher order terms, other transformations etc.
# nsp         - integer giving the number of starting points in each dimension to be 
#               used by optim. A lattice of points will be created by constructing 
#               an equally spaced grid for each of the (p - 1) hyperspherical coordinates
#               used to represent theta in the optimization. Default is 3
# L           - a list of integers specifying which starting points to use. If L = 0 (default),
#               all of the starting points in the lattice will be utilized. Otherwise,
#               L of these will be chosen by row number. If L = -1, the user will have to 
#               input a matrix whose rows are the starting points.
# etaStart    - a matrix with (p-1) columns each row indicating a unique starting value
#               used in optimization for estimating theta. This is input only if L=-1 
# datosx      - the dataset of n whose columns include the covariates, survey variables of the model.
# sp          - order of spline.
# dfs         - degrees of freedem of the spline

## Output: A List with the following elements
#
# thetaHat - length p vector giving the estimated coefficient
# fnvalue  - achieved minimum value of the criterion function for estimating theta
# etaStart - matrix with (p - 1) columns, each row indicating a unique starting value
#            used in optimization for estimating theta
# optInf   - list containing information about optimization routine for each
#            starting point

PLFSI_model <- function(si_vars = NULL, linear_vars = NULL, formula_lv=NULL, 
                  datosfda = NULL, tt = NULL, datosx=NULL, nsp = 3, L = 0, 
                  etaStart = NULL, sp=NULL, dfs=NULL) {
  
  library("numbers")
  # Perform checks
  if(is.null(si_vars)){
    stop('Must provide covariates for the single index part of the model')
  }
  if(is.null(linear_vars)){
    stop('Must provide covariate names for the linear part of the model')
  }
  if(is.null(datosx)){
    stop('Must provide covariate data for the model')
  }
  if(is.null(tt)){
    stop('Must provide grid vector tt for quantile functions')
  }
  if(is.null(datosfda) | !is.matrix(datosfda)){
    stop('Must provide responses as a matrix')
  }
  if(nrow(datosfda) != nrow(datosx) | ncol(datosfda) != length(tt)){
    stop('Dimensions of response do not match with tt and/or covariates')
  }
  if(!is.numeric(nsp) | length(nsp) != 1 | mod(nsp, 1) != 0 | nsp <= 0){
    message('Invalid specification of input nsp, resetting to default')
    nsp <- 3
  }
  if(!is.numeric(L) | length(L) != 1 | mod(L, 1) != 0){
    message('Invalid specification of input L, resetting to default')
    nsp <- L
  }
  if(L== -1 & is.null(etaStart)){
    message('Invalid specification of input etaStart, resetting to default')
    nsp <- L
  }
  
  # Create grid of starting values for optimization
  
  if(L!=-1) {
  # compute dimension of covariate in single index
  p <- length(si_vars)
  # specified spacing between staring points in each coordinate
  spc <- pi/nsp  
  # equally spaced starting points in polar coordinates
  f <- lapply(1:(p - 1), function(j) seq(-pi/2 + spc/2, pi/2 - spc/2, by = spc)) 
  # create grid of starting values
  etaStart <- as.matrix(expand.grid(f))   
  if(L != 0) {
    smp <- sample.int(nrow(etaStart), min(L, nrow(etaStart)))
    etaStart <- etaStart[smp,]
  }
  }
  ## To provide information about optimization as output
  optInf <- list()
  
  # provide criteria for termination of the algorithm
  optim_optns <- list(factr = 1e11, maxit = 100)
  
  WnMin <- rep(NA, nrow(as.matrix(etaStart)))
  etaMin <- matrix(NA, nrow = nrow(as.matrix(etaStart)), ncol = p - 1)
  converge<- rep(NA, nrow(as.matrix(etaStart)))
  
  #negll <- function(par, x, sleep=0, verbose=TRUE){
  #  if(verbose)
  #    cat(par, "\n")
  #  Sys.sleep(sleep)
  #  -sum(dnorm(x=x, mean=par[1], sd=par[2], log=TRUE))
  #}
  
  cl <- makeCluster(2)     # set the number of processor cores
  setDefaultCluster(cl=cl)
  
  source("wn_cost.R", local= knitr::knit_global())
  # main optimization loop over starting values
  for(k in 1:nrow(as.matrix(etaStart))) {
    print(k)
    WnOpt <- optimParallel(par = as.matrix(etaStart)[k,], fn = wn_cost, method = "L-BFGS-B",
                           lower = -pi/2, upper = pi/2, control = optim_optns, 
                           datosfda=datosfda, datosx = datosx, linear_vars=linear_vars, 
                           si_vars=si_vars, tt=tt, sp=sp, dfs=dfs, formula_lv=formula_lv,
                           parallel=list(forward=TRUE))
    
    optInf[[k]] <- WnOpt
    
    WnMin[k] <- WnOpt$value
    etaMin[k,] <- WnOpt$par
    converge[k]<- WnOpt$convergence
  }
  
  # the optimizer, i.e. thetaHat 
  source("polar2cart.R", local = knitr::knit_global())
  thetaHat <- polar2cart(etaMin[which.min(WnMin),], 1)
  
  optvalue <- min(WnMin)  # updated to find the minimized Wn in training set
  
  bs_fun <- as.matrix(datosx[,si_vars])%*%thetaHat
  
  bs_data<- data.frame(bs(x=bs_fun, degree=sp, df=dfs))
  colnames(bs_data) <- paste("BS", seq(1:ncol(bs_data)), sep = "")
  data_temp<- cbind.data.frame(datosx, bs_data)
  
  ## modifying the inputs of the algorithm
  
  formula<- paste(formula_lv, paste(colnames(bs_data),collapse="+"),sep="+")
  data_temp<- cbind.data.frame(data_temp, datosfda)
  data_analysis_svy <- svydesign(id= ~survey_id, strata = ~survey_strata, 
                                 weights= ~survey_wt, data = data_temp, nest = TRUE)
  
  objetofda= fdata(datosfda,argvals = tt)
  source("survey2wassersteinmodel.R", local =knitr::knit_global())
  res<- survey2wassersteinmodel(formula=formula, data_analysis_svy, objetofda=objetofda)
  
  q_n<- nrow(res$betaj)
  source("adj_fr_r2.R", local = knitr::knit_global())
  adjusted_fr_r2<- adj_fr_r2(qin=datosfda, qpred=res$predicciones, tt=tt, q=q_n, survey_weights= datosx$survey_wt)
  return(list(thetaHat = thetaHat, converge= converge, fnvalue=optvalue, etaStart= etaStart, optInf = optInf,
              R2= res$r2, betaj= res$betaj, predictions= res$predicciones, residuals= res$residuos,
              Frechet_R2= adjusted_fr_r2$Frechet_R2, Adj_Frechet_R2= adjusted_fr_r2$Adj_Frechet_R2,
              R2_vector= res$R2_vector))
  
}
