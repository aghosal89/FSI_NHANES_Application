

## This is the cost function for estimating the index parameter

## Inputs: 1) et           - index parameter in polar coordinates
##         2) datosfda     - an nxm matrix as the response whose each row represents an observation, 
##                           each column represents a quantile, t from [0,1].
##         3) categorical  - the names of categorical variables for linear part of model. 
##         4) numerical    - the names of numerical variables considered for single index model. 
##         5) tt           - the equidistant grid on [0,1] of length m.
##         6) sp           - the degree of polynomial considered for spline regression.
##         7) datos        - the dataset contaning all the covariates and the response. 
##         8) dfs          - degrees of freedom as an alternative to specifying the knots. 

## Output:the mean square prediction error. 




wn_cost <- function(et, datosfda, categorical, numerical, tt, sp, datos, dfs) {
  library("survey")
  library("splines")
  library("fda.usc")
  
  source("polar2cart.R", local= knitr::knit_global())
  th<- matrix(polar2cart(et, 1), length(et)+1, 1)
  bs_fun<- as.matrix(datos[,numerical])%*%th
  
  bs_data<- data.frame(bs(x=bs_fun, degree=sp, df=dfs))
  colnames(bs_data)<- paste("BS", seq(1:ncol(bs_data)), sep ="")
  data_temp<- data.frame(cbind.data.frame(datos, bs_data))
  
  data_analysis_svy <- svydesign(id= ~SDMVPSU, strata = ~SDMVSTRA, weights= ~wtmec4yr_adj_norm,
                                 data = data_temp,nest = TRUE)
  data_analysis_svy$variables$RIDRETH3= as.factor(data_analysis_svy$variables$RIDRETH3)
  data_analysis_svy$variables$RIAGENDR= as.factor(data_analysis_svy$variables$RIAGENDR)
  #data_analysis_svy$variables$edadcategorica= as.factor(data_analysis_svy$variables$edadcategorica)
  
  objetofda= fdata(datosfda,argvals = tt)
  formula<- paste(paste(categorical,collapse="+"),paste(colnames(bs_data),collapse="+"),sep="+")
  source("survey2wassersteinmodel.R", local =knitr::knit_global())
  res=  survey2wassersteinmodel(formulas=formula,data_analysis_svy, objetofda = objetofda)
  Yhat<- res$predicciones
  gx<- apply(Yhat, 2, mean)*dim(Yhat)[1]
  source("cuadratico.R", local =knitr::knit_global())
  if(any(diff(gx) < 0)) {
    Yhat<- cuadratico(Yhat)
  }
  
  #get array of Wasserstein distances between the response and its prediction
  RSS <- sapply(1:nrow(datosfda), function(i) fdadensity:::trapzRcpp(X = tt, Y = (datosfda[i, ] - Yhat[i, ])^2))
  
  w<- data_analysis_svy$variables$wtmec4yr_adj_norm
  return(sum(w*RSS/sum(w)))
}

